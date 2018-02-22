var ViewTree = (function() {
  var ViewTree = {};

  ViewTree.getView = function() {
    return layout.ViewTree;
  };

  ViewTree.initial = function() {
    ViewTree.initialized = true;
    var elem = ViewTree.getView();
    ViewTree.clock = new THREE.Clock();
    ViewTree.scene = new THREE.Scene();
    ViewTree.targetList = [];
    var width = $(elem).width();
    var height = $(elem).height();
    ViewTree.camera = new THREE.PerspectiveCamera(60, width / height, 0.1, 100);
    ViewTree.camera.position.z = 3;
    ViewTree.renderer = new THREE.WebGLRenderer({preserveDrawingBuffer: true});
    ViewTree.renderer.setPixelRatio(window.devicePixelRatio);
    ViewTree.renderer.setSize(width, height);
    ViewTree.renderer.setClearColor(0xffffff, 1);
    elem.appendChild(ViewTree.renderer.domElement);
    ViewTree.renderer.domElement.id = 'view-tree-canvas';

    var pointLight = new THREE.PointLight(0xffffff);
    pointLight.position.x = 100;
    pointLight.position.y = 100;
    pointLight.position.z = 100;
    ViewTree.scene.add(pointLight);

    ViewTree.directionalLight = new THREE.DirectionalLight(0xffffff);
    ViewTree.scene.add(ViewTree.directionalLight);

    var controls = new THREE.OrbitControls(ViewTree.camera, ViewTree.renderer.domElement);
    controls.mouseButtons = {
      ZOOM: THREE.MOUSE.MIDDLE,
      PAN: THREE.MOUSE.RIGHT,
      ORBIT: THREE.MOUSE.LEFT
    }
    ViewTree.raycaster = new THREE.Raycaster();

    ViewTree.render();
  };

  ViewTree.render = function() {
    // ViewTree.raycaster.setFromCamera(mousePos, camera);

    // scene
    var delta = ViewTree.clock.getDelta();
    requestAnimationFrame(ViewTree.render);
    // cameraControls.update(delta);
    ViewTree.directionalLight.position.copy(ViewTree.camera.position);
    ViewTree.renderer.render(ViewTree.scene, ViewTree.camera);
  };

  ViewTree.drawRadialTree = function() {
    var PARTICLE_SIZE = 20;
    var pointPositions = [];
    var pointColors = [];
    var linePositions = [];
    var lineColors = [];

    for (var id in data.treeData) {
      var i = Number(id);
      
      // pointPositions.push(data.treeData[id].x, data.treeData[id].y, 0);
      pointPositions.push(data.treeData[id].x, 0, data.treeData[id].y);
      pointColors.push( 0.01 + 0.1 * ( i / data.treeNodeCount ), 1.0, 0.5 );
      var parent = data.treeData[id].parent;
      if (parent != undefined) {
        var parentId = parent.id;
        // linePositions.push(data.treeData[id].x, data.treeData[id].y, 0);
        linePositions.push(data.treeData[id].x, 0, data.treeData[id].y);
        // linePositions.push(data.treeData[parentId].x, data.treeData[parentId].y, 0);
        linePositions.push(data.treeData[parentId].x, 0, data.treeData[parentId].y);
        lineColors.push(0.01 + 0.1 * ( i / data.treeNodeCount ), 1.0, 0.5);
        lineColors.push(0.01 + 0.1 * ( i / data.treeNodeCount ), 1.0, 0.5);
      }
    }
    console.log(pointPositions, pointColors);
    var pointGeometry = new THREE.BufferGeometry();
    pointGeometry.addAttribute('position', new THREE.Float32BufferAttribute(pointPositions, 3));
    pointGeometry.addAttribute('color', new THREE.Float32BufferAttribute(pointColors, 3));
    pointGeometry.computeBoundingSphere();
    var pointSprite = new THREE.TextureLoader().load("img/disc.png");
    var pointMaterial = new THREE.PointsMaterial({
      size: .01, 
      vertexColors: THREE.VertexColors, 
      map: pointSprite,
      alphaTest: 0.5,
      transparent: true
    });
    var points = new THREE.Points( pointGeometry, pointMaterial );
    console.log(points);
    ViewTree.scene.add( points );

    var lineGeometry = new THREE.BufferGeometry();
    var lineGeometry = new THREE.BufferGeometry();
    lineGeometry.addAttribute( 'position', new THREE.Float32BufferAttribute( linePositions, 3 ) );
    lineGeometry.addAttribute( 'color', new THREE.Float32BufferAttribute( lineColors, 3 ) );
    var lineMaterial = new THREE.LineBasicMaterial({color: 0x0000ff});
    var lines = new THREE.LineSegments(lineGeometry, lineMaterial);
    console.log(lines);
    ViewTree.scene.add( lines );
    
  }

  ViewTree.resize = function() {
    if (!ViewTree.initialized) return;
    var elem = ViewTree.getView();
    var width = $(elem).width();
    var height = $(elem).height();
    ViewTree.camera.aspect = width / height;
    ViewTree.camera.updateProjectionMatrix();
    ViewTree.renderer.setSize(width, height);
  }

  function getCameraMatrix() {
    return ViewTree.camera.matrix.toArray();
  }

  function setCameraMatrix(str) {
    var cameraState = JSON.parse(str);
    // ... read cameraState somehow ...
    ViewTree.camera.matrix.fromArray(cameraState);
    // Get back position/rotation/scale attributes
    ViewTree.camera.matrix.decompose(ViewTree.camera.position, ViewTree.camera.quaternion, ViewTree.camera.scale); 
  }

  return ViewTree;
})();
