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

var valueObject, labelObject;
function updateMesh(mesh) {
  data.mesh = mesh;
  console.log("updating mesh...");
  var valueGeometry = new THREE.Geometry();
  var labelGeometry = new THREE.Geometry();
  for (var i = 0; i < data.mesh.nNodes; i ++) {
    valueGeometry.vertices.push(new THREE.Vector3(data.mesh.coords[i*2] - 1.7, data.mesh.coords[i*2+1], 0));
    labelGeometry.vertices.push(new THREE.Vector3(data.mesh.coords[i*2] - 1.7, data.mesh.coords[i*2+1], 0));
  }
  for (var i = 0; i < data.mesh.nTriangles; i ++) {
    var face = new THREE.Face3(data.mesh.conn[i*3], data.mesh.conn[i*3+1], data.mesh.conn[i*3+2], null);
    valueGeometry.faces.push(face);
    labelGeometry.faces.push(face);
  }

  var valueMaterial = new THREE.MeshBasicMaterial( {vertexColors: THREE.FaceColors, side: THREE.DoubleSide} );
  var labelMaterial = new THREE.MeshBasicMaterial( {vertexColors: THREE.FaceColors, side: THREE.DoubleSide} );
  valueObject = new THREE.Mesh(valueGeometry, valueMaterial);
  labelObject = new THREE.Mesh(labelGeometry, labelMaterial);
  // scene.add(valueObject);
  scene.add(labelObject);
  requestData();
}

function updateData(values, labels) {
  data.values = values;
  data.labels = labels;
  var range = {};
  range.max = Math.max(...data.values);
  range.min = Math.min(...data.values);
  
  console.log(range);
  var count1 = 0;
  var count2 = 0;
  for (var i = 0; i < data.values.length; i ++) {
    if (data.values[i] < 0) count2 ++;
    if (data.values[i] > 0) count1 ++;
  }
  console.log(count1, count2);

  var cs = d3.scaleLinear().domain([range.min, 0, range.max])
      .range([d3.rgb("#ff0000"), d3.rgb('#ffffff'), d3.rgb('#0000ff')]);
  var valueColorScale = function(v) {
    return new THREE.Color(cs(v).toString());
  };
  for (var i = 0; i < valueObject.geometry.faces.length; i ++) {

    // update value color
    var face = valueObject.geometry.faces[i];
    var v1 = data.values[face.a];
    var v2 = data.values[face.b];
    var v3 = data.values[face.c];
    if (valueObject.geometry.faces[i].vertexColors.length == 0) {
      valueObject.geometry.faces[i].vertexColors.push(valueColorScale(v1));
      valueObject.geometry.faces[i].vertexColors.push(valueColorScale(v2));
      valueObject.geometry.faces[i].vertexColors.push(valueColorScale(v3));
    }
    else {
      valueObject.geometry.faces[i].vertexColors[0].set(valueColorScale(v1));
      valueObject.geometry.faces[i].vertexColors[1].set(valueColorScale(v2));
      valueObject.geometry.faces[i].vertexColors[2].set(valueColorScale(v3));
    }

    // update label color
    var l1 = data.labels[face.a];
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
    labelObject.geometry.faces[i].color.set(labelColor);
  }
  valueObject.geometry.colorsNeedUpdate = true;
  valueObject.material.needsUpdate = true;
  labelObject.geometry.colorsNeedUpdate = true;
  labelObject.material.needsUpdate = true;
  $('#loading').hide();
}

function updateRenderMethod(method) {
  console.log(method);
  if (method === 'label') {
    scene.remove(valueObject);
    scene.add(labelObject);
  }
  else if (method === 'value'){
    scene.remove(labelObject);
    scene.add(valueObject);
  }
}

function updateRenderWireframe(renderWireframe) {
  valueObject.material.wireframe = renderWireframe;
  valueObject.geometry.colorsNeedUpdate = true;
  valueObject.material.needsUpdate = true;
}

initializeControlPanel();
render();
