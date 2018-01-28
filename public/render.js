var stats = new Stats();
stats.showPanel(0);
stats.dom.id="stats";
document.body.appendChild(stats.dom);
$("#stats").css({visibility: "hidden"});

var clock = new THREE.Clock();
var scene = new THREE.Scene();

var camera = new THREE.PerspectiveCamera(30, window.innerWidth/window.innerHeight, 0.1, 100); 
camera.position.z = 3;

var renderer = new THREE.WebGLRenderer({preserveDrawingBuffer: true});
renderer.setPixelRatio(window.devicePixelRatio);
renderer.setSize(window.innerWidth, window.innerHeight);
renderer.setClearColor(0xffffff, 1);
document.body.appendChild(renderer.domElement);

var pointLight = new THREE.PointLight(0xffffff);
pointLight.position.x = 100;
pointLight.position.y = 100;
pointLight.position.z = 100;
// pointLight.castShadow = true;
// pointLight.shadowDarkness = 0.5;
scene.add(pointLight);

var directionalLight = new THREE.DirectionalLight(0xffffff);
scene.add(directionalLight);

cameraControls = new THREE.TrackballControls(camera, renderer.domElement);
cameraControls.target.set(0, 0, 0);
cameraControls.zoomSpeed = 0.04;
cameraControls.panSpeed = 0.04;
// cameraControls.addEventListener("change", render); // not working.. sigh

var raycaster = new THREE.Raycaster();
var mousePos = new THREE.Vector2();

window.addEventListener("mousedown", onMouseDown, false );
window.addEventListener("mousemove", onMouseMove, false);
window.addEventListener("resize", onResize, false);

function render() {
  stats.begin();

  raycaster.setFromCamera(mousePos, camera);
  var intersects = raycaster.intersectObjects(scene.children);
  // for (i=0; i<intersects.length; i++)
  //   intersects[i].object.material.color.set(0xff0000);

  // scene
  var delta = clock.getDelta();
  requestAnimationFrame(render);
  cameraControls.update(delta);
  directionalLight.position.copy(camera.position);
  renderer.render(scene, camera);

  stats.end();
}

function onMouseDown(evt) {
  mousePos.x = (evt.clientX / window.innerWidth) * 2 - 1;
  mousePos.y = -(evt.clientY / window.innerHeight) * 2 + 1;
}

function onMouseMove(evt) {
  mousePos.x = (evt.clientX / window.innerWidth) * 2 - 1;
  mousePos.y = -(evt.clientY / window.innerHeight) * 2 + 1;
}

function onResize() {
  camera.aspect = window.innerWidth/window.innerHeight;
  camera.updateProjectionMatrix();
  renderer.setSize(window.innerWidth, window.innerHeight);
  cameraControls.handleResize();
}

var geom;
var material;
function updateMesh(mesh) {
  data.mesh = mesh;
  console.log("updating mesh...");
  geom = new THREE.Geometry();
  for (var i = 0; i < data.mesh.nNodes; i ++) {
    geom.vertices.push(new THREE.Vector3(data.mesh.coords[i*2] - 1.7, data.mesh.coords[i*2+1], 0));
  }
  for (var i = 0; i < data.mesh.nTriangles; i ++) {
    var face = new THREE.Face3(data.mesh.conn[i*3], data.mesh.conn[i*3+1], data.mesh.conn[i*3+2], null);// , new THREE.Color(Math.random() * 0xfffffff));
    geom.faces.push(face);
    face.color.set(new THREE.Color(Math.random() * 0xfffffff));
  }

  material = new THREE.MeshBasicMaterial( {vertexColors: THREE.VertexColors, side: THREE.DoubleSide} );
  var obj = new THREE.Mesh( geom, material );
  scene.add(obj);

  requestData();
}

function updateData(values, labels) {
  data.values = values;
  data.labels = labels;
  var range = {};
  range.max = Math.max(...data.values);
  range.min = Math.min(...data.values);
  var colorScale = function(v, positive, negative, zero) {
    var c = zero;
    if (v > 0) c = positive * v / range.max;
    if (v < 0) c = negative * v / range.min;
    return new THREE.Color(c);
  };
  for (var i = 0; i < geom.faces.length; i ++) {
    var face = geom.faces[i];
    /* var l1 = data.labels[face.a];
    var l2 = data.labels[face.b];
    var l3 = data.labels[face.c];
    var labelColor = new THREE.Color(0x000000);
    if (l1 === l2 && l2 === l3) {
      if (l1 === 0) {
        labelColor = new THREE.Color(0x000000);
      }
      else if (l1 === 1) {
        labelColor = new THREE.Color(0xff0000);
      }
      else {
        labelColor = new THREE.Color(0x0000ff);
      }
    }
    geom.faces[i].color.set(labelColor); */

    var v1 = data.values[face.a];
    var v2 = data.values[face.b];
    var v3 = data.values[face.c];
    var positive = 0x00ff00;
    var negative = 0xff0000;
    var zero = 0xffffff;
    geom.faces[i].vertexColors[0] = colorScale(v1, positive, negative, zero);
    geom.faces[i].vertexColors[1] = colorScale(v2, positive, negative, zero);
    geom.faces[i].vertexColors[2] = colorScale(v3, positive, negative, zero);
  }
  geom.colorsNeedUpdate = true;
  material.vertexColors = THREE.VertexColors;
  material.needsUpdate = true;
}

initializeControlPanel();
render();
