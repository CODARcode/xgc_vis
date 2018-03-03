var View3D = (function() {
  var View3D = {};

  View3D.getView = function() {
    return layout.View3D;
  };

  View3D.getWH = function() {
    return {
      width: layout.View3D.clientWidth, 
      height: layout.View3D.clientHeight
    };
  };

  View3D.getMatrix = function() {
    var projectionMatrix = View3D.camera.projectionMatrix.clone();
    var cameraModelMatrix = View3D.camera.matrixWorld.clone();
    var cameraModelMatrixInverse = View3D.camera.matrixWorldInverse.clone();
    var mul = projectionMatrix.multiply(cameraModelMatrixInverse);
    var inv = new THREE.Matrix4();
    inv.getInverse(mul);
    return inv.elements;
  };

  View3D.getCameraFov = function() {
    return View3D.camera.fov;
  };

  View3D.getCameraAspectRatio = function() {
    return View3D.camera.aspect;
  };

  View3D.getCameraNear = function() {
    return View3D.camera.near;
  };

  View3D.getCameraFar = function() {
    return View3D.camera.far;
  };

  View3D.getCameraMatrix = function() {
    return View3D.camera.matrix.toArray();
  };

  View3D.setCameraMatrix = function(str) {
    var cameraState = JSON.parse(str);
    View3D.camera.matrix.fromArray(cameraState);
    View3D.camera.matrix.decompose(View3D.camera.position, View3D.camera.quaternion, View3D.camera.scale); 
  };

  View3D.resetCamera = function() {
    View3D.controls.reset();
  };

  View3D.getRenderer = function() {
    return View3D.renderer;
  };

  View3D.updateLegendDisplay = function(display) {
    if (display) {
      $('.legend-div-3d').show();
    }
    else {
      $('.legend-div-3d').hide();
    }
  };

  View3D.updateTimestepDisplay = function(display) {
    if (display) {
      $('.timestep-div-3d').show();
    }
    else {
      $('.timestep-div-3d').hide();
    }
  };

  View3D.updateStartAngle = function(startAngle) {

  };

  View3D.initial = function() {
    View3D.initialized = true;
    var elem = View3D.getView();
    View3D.clock = new THREE.Clock();
    View3D.scene = new THREE.Scene();
    View3D.targetList = [];
    var width = $(elem).width();
    var height = $(elem).height();
    // View3D.camera = new THREE.PerspectiveCamera(30, width / height, 0.1, 100);
    View3D.camera = new THREE.PerspectiveCamera(30, width / height, 0.1, 100);
    // View3D.camera.position.z = 3;
    View3D.camera.position.set(0, 0, 6);
    View3D.renderer = new THREE.WebGLRenderer({preserveDrawingBuffer: true});
    View3D.renderer.setPixelRatio(window.devicePixelRatio);
    View3D.renderer.setSize(width, height);
    View3D.renderer.setClearColor(0xffffff, 1);
    elem.appendChild(View3D.renderer.domElement);
    View3D.renderer.domElement.id = 'view-3d-canvas';

    var pointLight = new THREE.PointLight(0xffffff);
    pointLight.position.x = 100;
    pointLight.position.y = 100;
    pointLight.position.z = 100;
    View3D.scene.add(pointLight);

    View3D.directionalLight = new THREE.DirectionalLight(0xffffff);
    View3D.scene.add(View3D.directionalLight);

    // View3D.controls = new THREE.OrbitControls(View3D.camera, View3D.renderer.domElement);
    // View3D.controls.mouseButtons = {
    //   ZOOM: THREE.MOUSE.MIDDLE,
    //   PAN: THREE.MOUSE.RIGHT,
    //   ORBIT: THREE.MOUSE.LEFT
    // }
    // View3D.controls.minAzimuthAngle = - Infinity;
    // View3D.controls.maxAzimuthAngle = Infinity;
    // View3D.controls.maxPolarAngle = Infinity;
    // View3D.controls.minPolarAngle = -Infinity;

    View3D.controls = new THREE.TrackballControls(View3D.camera, View3D.getView());
    View3D.controls.rotateSpeed = 1.1;
    View3D.controls.zoomSpeed = .1;
    View3D.controls.panSpeed = 0.8;
    View3D.controls.noZoom = false;
    View3D.controls.noPan = false;
    View3D.controls.staticMoving = true;
    View3D.controls.dynamicDampingFactor = 0.3;
    View3D.controls.keys = [65, 83, 68];
    View3D.controls.addEventListener('change', View3D.render);

    View3D.raycaster = new THREE.Raycaster();

    View3D.render();
    View3D.animate();

    var imgNode = $('.volren-image').detach();
    var parentNode = View3D.getView();
    $(parentNode).append(imgNode);

    $(parentNode).mouseup(function() {
      Client.requestImage();
    });

    $(window).bind('mousewheel', function(e){
      var elem = View3D.getView();
      if (e.clientX > elem.offsetLeft && e.clientX < elem.offsetLeft + elem.offsetWidth
        && e.clientY > elem.offsetTop && e.clientY < elem.offsetTop + elem.offsetHeight) {
        View3D.hideVolrenImage();
        Client.requestImageWait();
      }
    });

    $(parentNode).bind('touchmove', function(e) {
      View3D.hideVolrenImage();
      Client.requestImageWait();
    });

    $(parentNode).mousedown(function() {
      View3D.hideVolrenImage();
    });
  };

  View3D.hideVolrenImage = function() {
    $('.volren-image').hide();
  };

  View3D.showVolrenImage = function() {
    $('.volren-image').show();
  };

  View3D.updateImage = function(imageBlob) {
    var objectURL = URL.createObjectURL(imageBlob);
    var img = document.querySelector('.volren-image');
    img.src = objectURL;
    var imgNode = $('.volren-image')[0];
    var parentNode = View3D.getView();
    var top = '-' + $(parentNode).height() + 'px';
    $(imgNode).css('top', top);
    View3D.showVolrenImage();
    // $('.volren-image').html('<img src="data:image/png;base64,' + data + '" />');
  };

  View3D.drawMesh = function() {
    while(View3D.scene.children.length > 0){ 
      View3D.scene.remove(View3D.scene.children[0]); 
    }
    console.debug("draw 3D meshes...");

    View3D.valueObjectList = [];
    View3D.labelObjectList = [];
    for (var phiIndex = 0; phiIndex < data.nPhi; phiIndex ++) {
      var valueGeometry = new THREE.Geometry();
      var labelGeometry = new THREE.Geometry();
      for (var i = 0; i < data.mesh.nNodes; i ++) {
        labelGeometry.vertices.push(new THREE.Vector3(data.mesh.coords[i*2], 0, data.mesh.coords[i*2+1]));
        valueGeometry.vertices.push(new THREE.Vector3(data.mesh.coords[i*2], 0, data.mesh.coords[i*2+1]));
      }
      for (var i = 0; i < data.mesh.nTriangles; i ++) {
        labelGeometry.faces.push(new THREE.Face3(data.mesh.conn[i*3], data.mesh.conn[i*3+1], data.mesh.conn[i*3+2], null));
        valueGeometry.faces.push(new THREE.Face3(data.mesh.conn[i*3], data.mesh.conn[i*3+1], data.mesh.conn[i*3+2], null));
        valueGeometry.faces[i].vertexColors.push(new THREE.Color(0x00ffff));
        valueGeometry.faces[i].vertexColors.push(new THREE.Color(0x00ffff));
        valueGeometry.faces[i].vertexColors.push(new THREE.Color(0x00ffff));
      }

      var labelMaterial = new THREE.MeshBasicMaterial( {vertexColors: THREE.FaceColors, side: THREE.DoubleSide} );
      var valueMaterial = new THREE.MeshBasicMaterial( {vertexColors: THREE.FaceColors, side: THREE.DoubleSide} );
      var labelObject = new THREE.Mesh(labelGeometry, labelMaterial);
      var valueObject = new THREE.Mesh(valueGeometry, valueMaterial);

      var angle = Math.PI * 2 * phiIndex / data.nPhi;
      labelObject.rotateZ(angle);
      valueObject.rotateZ(angle);
      View3D.valueObjectList.push(valueObject);
      View3D.labelObjectList.push(labelObject);

      if (!menuText.renderMethod || menuText.renderMethod === 'value') {
        View3D.targetList.concat(valueObject);
        View3D.scene.add(valueObject);
      }
      else {
        View3D.targetList.concat(labelObject);
        View3D.scene.add(labelObject);
      }
    }
  };

  View3D.render = function() {
    var delta = View3D.clock.getDelta();
    View3D.directionalLight.position.copy(View3D.camera.position);
    View3D.renderer.render(View3D.scene, View3D.camera);
  };

  View3D.animate = function() {
    requestAnimationFrame(View3D.animate);
    View3D.controls.update();
    View3D.render();
  };

  View3D.resize = function() {
    if (!View3D.initialized) return;
    var elem = View3D.getView();
    var width = $(elem).width();
    var height = $(elem).height();
    View3D.camera.aspect = width / height;
    View3D.camera.updateProjectionMatrix();
    View3D.renderer.setSize(width, height);
    View3D.controls.handleResize();
    ViewTree.render();
    View3D.hideVolrenImage();
    if (ViewTF.initialized) {
      Client.requestImageWait();
    }
  }

  View3D.updateData = function() {
    console.debug('drawing all slices data');
    console.log(data.slices);
    var red = '#b82e2e';
    var green = '#109618';
    var blue = '#3366cc';
    var white = '#ffffff';
    var black = '#000000';
    var redHex = 0xb82e2e;
    var greenHex = 0x109618;
    var blueHex = 0x3366cc;
    var whiteHex = 0xffffff;
    var blackHex = 0x000000;
    
    /* value color scale */
    var range = { max : 0, min : 0};
    data.slices.forEach(function(arr) {
      arr.forEach(function(d) {
        if (d > range.max) range.max = d;
        if (d < range.min) range.min = d;
      });
    });
    console.log(range);
    var cs = d3.scaleLinear().domain([range.min, 0, range.max])
        .range([d3.rgb(blue), d3.rgb(white), d3.rgb(red)]);
    var alphaScaleN = d3.scaleLinear().domain([range.min, 0]).range(0, 1);
    var alphaScaleP = d3.scaleLinear().domain([0, range.max]).range(0, 1);
    var valueColorScale = function(v) {
      var color = new THREE.Color(cs(v).toString());
      var alpha = v > 0 ? alphaScaleP(v) : alphaScaleN(v);
      color.set
      return color;
    };

    /* value color legend */
    var formatToOne = d3.format(".1e");
    d3.select('#min-value-text').text(formatToOne(range.min));
    d3.select('#max-value-text').text(formatToOne(range.max));

    /* label color scale */
    var labelRange = {};
    data.slices.forEach(function(arr) {
      arr.forEach(function(d) {
        if (d > labelRange.max) labelRange.max = d;
        if (d < labelRange.min) labelRange.min = d;
      });
    });
    // console.log(labelRange);
    var cs2 = d3.scaleLinear().domain([labelRange.min, 0, labelRange.max])
        .range([d3.rgb(green), d3.rgb(black), d3.rgb(red)]);
    var labelColorScale = function(l) {
      return new THREE.Color(cs2(l).toString());
    };

    for (var phiIndex = 0; phiIndex < View3D.valueObjectList.length; phiIndex ++) {
      for (var i = 0; i < View3D.valueObjectList[phiIndex].geometry.faces.length; i ++) {
        var face = View3D.valueObjectList[phiIndex].geometry.faces[i];
        // update label color
        // var l1 = data.labels[face.a];
        // var l2 = data.labels[face.b];
        // var l3 = data.labels[face.c];
        // var labelColor = new THREE.Color(blackHex);
        // if (l1 === l2 && l2 === l3) {
        //   // labelColor = labelColorScale(l1);
        //   if (l1 === 0) {
        //     labelColor = new THREE.Color(blackHex);
        //   }
        //   else if (l1 > 0) {
        //     labelColor = new THREE.Color(redHex);
        //   }
        //   else {
        //     labelColor = new THREE.Color(blueHex);
        //   }
        // }
        // labelObject.geometry.faces[i].color.set(labelColor);

        // update value color
        var v1 = data.slices[phiIndex][face.a];
        var v2 = data.slices[phiIndex][face.b];
        var v3 = data.slices[phiIndex][face.c];
        View3D.valueObjectList[phiIndex].geometry.faces[i].vertexColors[0].set(valueColorScale(v1));
        View3D.valueObjectList[phiIndex].geometry.faces[i].vertexColors[1].set(valueColorScale(v2));
        View3D.valueObjectList[phiIndex].geometry.faces[i].vertexColors[2].set(valueColorScale(v3));
      }
      View3D.valueObjectList[phiIndex].geometry.colorsNeedUpdate = true;
      View3D.valueObjectList[phiIndex].material.needsUpdate = true;
      // labelObject.geometry.colorsNeedUpdate = true;
      // labelObject.material.needsUpdate = true;
    }
    
  };

  return View3D;
})();
