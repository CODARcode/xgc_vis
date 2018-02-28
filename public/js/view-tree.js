var ViewTree = (function() {
  var ViewTree = {};

  ViewTree.getView = function() {
    return layout.ViewTree;
  };

  ViewTree.getCameraMatrix = function() {
    return ViewTree.camera.matrix.toArray();
  }

  ViewTree.setCameraMatrix = function(str) {
    var cameraState = JSON.parse(str);
    ViewTree.camera.matrix.fromArray(cameraState);
    ViewTree.camera.matrix.decompose(ViewTree.camera.position, ViewTree.camera.quaternion, ViewTree.camera.scale); 
  }

  ViewTree.resetCamera = function() {
    ViewTree.controls.reset();
  }

  ViewTree.getRenderer = function() {
    return ViewTree.renderer;
  }

  ViewTree.initial = function() {
    ViewTree.initialized = true;
    var elem = ViewTree.getView();
    ViewTree.clock = new THREE.Clock();
    ViewTree.scene = new THREE.Scene();
    ViewTree.targetList = [];
    var width = $(elem).width();
    var height = $(elem).height();
    ViewTree.camera = new THREE.PerspectiveCamera(60, width / height, 0.1, 1000);
    ViewTree.camera.position.set(0, 0, 800);
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

    // ViewTree.controls = new THREE.OrbitControls(ViewTree.camera, ViewTree.renderer.domElement);
    // ViewTree.controls.mouseButtons = {
    //   ZOOM: THREE.MOUSE.MIDDLE,
    //   PAN: THREE.MOUSE.RIGHT,
    //   ORBIT: THREE.MOUSE.LEFT
    // }

    ViewTree.controls = new THREE.TrackballControls(ViewTree.camera, ViewTree.getView());
    ViewTree.controls.rotateSpeed = 1.0;
    ViewTree.controls.zoomSpeed = .1;
    ViewTree.controls.panSpeed = 0.8;
    ViewTree.controls.noZoom = false;
    ViewTree.controls.noPan = false;
    ViewTree.controls.staticMoving = true;
    ViewTree.controls.dynamicDampingFactor = 0.3;
    ViewTree.controls.keys = [65, 83, 68];
    ViewTree.controls.addEventListener('change', ViewTree.render);

    ViewTree.raycaster = new THREE.Raycaster();

    ViewTree.render();
    ViewTree.animate();
  };

  ViewTree.render = function() {
    var delta = ViewTree.clock.getDelta();
    ViewTree.directionalLight.position.copy(ViewTree.camera.position);
    ViewTree.renderer.render(ViewTree.scene, ViewTree.camera);
  };

  ViewTree.animate = function() {
    requestAnimationFrame(ViewTree.animate);
    ViewTree.controls.update();
    ViewTree.render();
  };

  ViewTree.drawRadialTree = function() {
    var PARTICLE_SIZE = 20;
    var pointPositions = [];
    var pointColors = [];
    var linePositions = [];
    var lineColors = [];

    var cs = d3.scaleLinear().domain([data.treeMinMax.extremumMin, 0, data.treeMinMax.extremumMax])
        .range([d3.rgb(PresetColor.blue), d3.rgb(PresetColor.black), d3.rgb(PresetColor.red)]);
    var valueColorScale = function(v) {
      return new THREE.Color(cs(v).toString());
    };

    for (var id in data.treeData) {
      var i = Number(id);
      
      // point 
      pointPositions.push(data.treeData[id].x, data.treeData[id].extremumY, data.treeData[id].y);
      // pointColors.push((data.treeData[id].depth / data.treeDepth), 1.0, 1.0);
      var extremumColorArr = valueColorScale(data.treeData[id].data.extremum_val).toArray();
      pointColors.push(extremumColorArr[0], extremumColorArr[1], extremumColorArr[2]);

      // point-parent line 
      var parent = data.treeData[id].parent;
      var lineColor = new THREE.Color(PresetColor.gray).toArray();
      if (parent != undefined) {
        var parentId = parent.id;
        
        var saddleColorArr = valueColorScale(data.treeData[id].data.saddle_val).toArray();
        var parentExtremumColorArr = valueColorScale(data.treeData[parentId].data.extremum_val).toArray();

        linePositions.push(data.treeData[id].x, data.treeData[id].extremumY, data.treeData[id].y);
        linePositions.push(data.treeData[id].x, data.treeData[id].saddleY, data.treeData[id].y);
        // lineColors.push(extremumColorArr[0], extremumColorArr[1], extremumColorArr[2]);
        // lineColors.push(saddleColorArr[0], saddleColorArr[1], saddleColorArr[2]);
        lineColors.push(lineColor[0], lineColor[1], lineColor[2]);
        lineColors.push(lineColor[0], lineColor[1], lineColor[2]);

        linePositions.push(data.treeData[id].x, data.treeData[id].saddleY, data.treeData[id].y);
        linePositions.push(data.treeData[parentId].x, data.treeData[id].saddleY, data.treeData[parentId].y);
        // lineColors.push(saddleColorArr[0], saddleColorArr[1], saddleColorArr[2]);
        // lineColors.push(saddleColorArr[0], saddleColorArr[1], saddleColorArr[2]);
        lineColors.push(lineColor[0], lineColor[1], lineColor[2]);
        lineColors.push(lineColor[0], lineColor[1], lineColor[2]);

        linePositions.push(data.treeData[parentId].x, data.treeData[id].saddleY, data.treeData[parentId].y);
        linePositions.push(data.treeData[parentId].x, data.treeData[parentId].extremumY, data.treeData[parentId].y);
        // lineColors.push(saddleColorArr[0], saddleColorArr[1], saddleColorArr[2]);
        // lineColors.push(parentExtremumColorArr[0], parentExtremumColorArr[1], parentExtremumColorArr[2]);
        lineColors.push(lineColor[0], lineColor[1], lineColor[2]);
        lineColors.push(lineColor[0], lineColor[1], lineColor[2]);
      }
    }
    console.log(pointPositions, pointColors);
    var pointGeometry = new THREE.BufferGeometry();
    pointGeometry.addAttribute('position', new THREE.Float32BufferAttribute(pointPositions, 3));
    pointGeometry.addAttribute('color', new THREE.Float32BufferAttribute(pointColors, 3));
    pointGeometry.computeBoundingSphere();
    var pointSprite = new THREE.TextureLoader().load("img/disc.png");
    var pointMaterial = new THREE.PointsMaterial({
      size: 10, 
      vertexColors: THREE.VertexColors, 
      map: pointSprite,
      alphaTest: 0.5,
      transparent: true
    });
    var points = new THREE.Points( pointGeometry, pointMaterial );
    console.log(points);
    ViewTree.scene.add( points );

    console.log(linePositions, lineColors);
    var lineGeometry = new THREE.BufferGeometry();
    lineGeometry.addAttribute('position', new THREE.Float32BufferAttribute(linePositions, 3));
    lineGeometry.addAttribute('color', new THREE.Float32BufferAttribute(lineColors, 3));
    var lineMaterial = new THREE.LineBasicMaterial({
      vertexColors: THREE.VertexColors, 
      opacity: 1
    });
    var lines = new THREE.LineSegments(lineGeometry, lineMaterial);
    console.log(lines);
    ViewTree.scene.add( lines );

    var zeroPlaneGeometry = new THREE.PlaneGeometry(50000, 50000);
    var zeroPlaneMaterial = new THREE.MeshBasicMaterial({color: 0xeeeeee, side: THREE.DoubleSide, transparent: true, opacity: .5});
    var zeroPlane = new THREE.Mesh(zeroPlaneGeometry, zeroPlaneMaterial);
    zeroPlane.rotateX(Math.PI / 2);
    ViewTree.scene.add(zeroPlane);
  }

  ViewTree.resize = function() {
    if (!ViewTree.initialized) return;
    var elem = ViewTree.getView();
    var width = $(elem).width();
    var height = $(elem).height();
    ViewTree.camera.aspect = width / height;
    ViewTree.camera.updateProjectionMatrix();
    ViewTree.renderer.setSize(width, height);
    ViewTree.controls.handleResize();
    ViewTree.render();
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
