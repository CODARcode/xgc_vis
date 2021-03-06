var View2D = (function() {
  var View2D = {};

  View2D.getView = function() {
    return layout.View2D;
  };

  View2D.getCameraMatrix = function() {
    return View2D.camera.matrix.toArray();
  };

  View2D.setCameraMatrix = function(str) {
    var cameraState = JSON.parse(str);
    View2D.camera.matrix.fromArray(cameraState);
    View2D.camera.matrix.decompose(View2D.camera.position, View2D.camera.quaternion, View2D.camera.scale); 
  };

  View2D.resetCamera = function() {
    View2D.controls.reset();
  };

  View2D.getRenderer = function() {
    return View2D.renderer;
  };

  View2D.updateLegendDisplay = function(display) {
    if (display) {
      $('.legend-div-2d').show();
    }
    else {
      $('.legend-div-2d').hide();
    }
  };

  View2D.updateTimestepDisplay = function(display) {
    if (display) {
      $('.timestep-div-2d').show();
    }
    else {
      $('.timestep-div-2d').hide();
    }
  };
  
  View2D.initial = function() {
    View2D.initialized = true;
    var elem = View2D.getView();
    View2D.clock = new THREE.Clock();
    View2D.scene = new THREE.Scene();
    View2D.targetList = [];
    var width = $(elem).width();
    var height = $(elem).height();
    // View2D.camera = new THREE.PerspectiveCamera(30, window.innerWidth/window.innerHeight, 0.1, 100);
    View2D.camera = new THREE.PerspectiveCamera(30, width / height, 0.1, 100);
    View2D.camera.position.z = 3;
    View2D.renderer = new THREE.WebGLRenderer({preserveDrawingBuffer: true});
    View2D.renderer.setPixelRatio(window.devicePixelRatio);
    // View2D.renderer.setSize(window.innerWidth, window.innerHeight);
    View2D.renderer.setSize(width, height);
    View2D.renderer.setClearColor(0xffffff, 1);
    elem.appendChild(View2D.renderer.domElement);
    View2D.renderer.domElement.id = 'view-2d-canvas';

    var pointLight = new THREE.PointLight(0xffffff);
    pointLight.position.x = 100;
    pointLight.position.y = 100;
    pointLight.position.z = 100;
    // pointLight.castShadow = true;
    // pointLight.shadowDarkness = 0.5;
    View2D.scene.add(pointLight);

    View2D.directionalLight = new THREE.DirectionalLight(0xffffff);
    View2D.scene.add(View2D.directionalLight);

    View2D.controls = new THREE.OrbitControls(View2D.camera, View2D.renderer.domElement);
    View2D.controls.enableRotate = false;
    View2D.controls.mouseButtons = {
      ZOOM: THREE.MOUSE.MIDDLE,
      PAN: THREE.MOUSE.LEFT
    }
    View2D.raycaster = new THREE.Raycaster();

    // window.addEventListener("resize", onResize, false);
    $(elem).mousemove(function(event) {
      onDocumentMouseMove(event);
    });

    $(elem).mouseout(function(evt) {
      $('.info-tooltip').hide();
    });

    View2D.render();
  };

  View2D.drawMesh = function() {
    while(View2D.scene.children.length > 0){ 
      View2D.scene.remove(View2D.scene.children[0]); 
    }
    console.debug("draw 2D mesh...");
    var valueGeometry = new THREE.Geometry();
    var labelGeometry = new THREE.Geometry();
    for (var i = 0; i < data.mesh.nNodes; i ++) {
      labelGeometry.vertices.push(new THREE.Vector3(data.mesh.coords[i*2] - 1.7, data.mesh.coords[i*2+1], 0));
      valueGeometry.vertices.push(new THREE.Vector3(data.mesh.coords[i*2] - 1.7, data.mesh.coords[i*2+1], 0));
    }
    for (var i = 0; i < data.mesh.nTriangles; i ++) {
      labelGeometry.faces.push(new THREE.Face3(data.mesh.conn[i*3], data.mesh.conn[i*3+1], data.mesh.conn[i*3+2], null));
      valueGeometry.faces.push(new THREE.Face3(data.mesh.conn[i*3], data.mesh.conn[i*3+1], data.mesh.conn[i*3+2], null));
      valueGeometry.faces[i].vertexColors.push(new THREE.Color(0xffffff));
      valueGeometry.faces[i].vertexColors.push(new THREE.Color(0xffffff));
      valueGeometry.faces[i].vertexColors.push(new THREE.Color(0xffffff));
    }

    var labelMaterial = new THREE.MeshBasicMaterial( {vertexColors: THREE.FaceColors, side: THREE.DoubleSide} );
    var valueMaterial = new THREE.MeshBasicMaterial( {vertexColors: THREE.FaceColors, side: THREE.DoubleSide} );
    View2D.labelObject = new THREE.Mesh(labelGeometry, labelMaterial);
    View2D.valueObject = new THREE.Mesh(valueGeometry, valueMaterial);
    if (!menuText.renderMethod || menuText.renderMethod === 'value') {
      View2D.targetList.push(View2D.valueObject);
      View2D.scene.add(View2D.valueObject);
    }
    else {
      View2D.targetList.push(View2D.labelObject);
      View2D.scene.add(View2D.labelObject);
    }

    var circleGeometry = new THREE.CircleGeometry(.001, 32);
    var circleMaterial = new THREE.MeshBasicMaterial( { color: 0xffff00 } );
    View2D.circleObject = new THREE.Mesh(circleGeometry, circleMaterial);
    View2D.scene.add(View2D.circleObject);
  };

  View2D.render = function() {
    // stats.begin();

    // View2D.raycaster.setFromCamera(View2D.mousePos, View2D.camera);

    // scene
    var delta = View2D.clock.getDelta();
    requestAnimationFrame(View2D.render);
    View2D.directionalLight.position.copy(View2D.camera.position);
    View2D.renderer.render(View2D.scene, View2D.camera);

    // stats.end();
  };

  View2D.resize = function() {
    if (!View2D.initialized) return;
    var elem = View2D.getView();
    var width = $(elem).width();
    var height = $(elem).height();
    // View2D.camera.aspect = window.innerWidth/window.innerHeight;
    View2D.camera.aspect = width / height;
    View2D.camera.updateProjectionMatrix();
    // View2D.renderer.setSize(window.innerWidth, window.innerHeight);
    View2D.renderer.setSize(width, height);
  }

  View2D.updateEnableSameTFEditor = function(enableSame) {
    globalStatus.updateEnableSameTFEditor = enableSame;
    if (enableSame) {
      View2D.updateTF(globalStatus.tfControlPoints);
    }
    else {
      View2D.drawData();
    }
  };

  View2D.updateTF = function(controlPoints) {
    if (!globalStatus.updateEnableSameTFEditor) {
      return;
    }
    if (!View2D.initialized || !View2D.dataInitialized) {
      return;
    }
    var domainList = [];
    var rangeList = [];
    controlPoints.forEach(function(d) {
      var v = data.range.min + (data.range.max - data.range.min) * d.alpha;
      v = Math.min(v, data.range.max);
      domainList.push(v);
      rangeList.push(d3.rgb(d.rgb));
    });
    var cs = d3.scaleLinear().domain(domainList)
        .range(rangeList);
    var valueColorScale = function(v) {
      return new THREE.Color(cs(v).toString());
    };
    View2D.drawDataColor(valueColorScale);
  };

  View2D.updateData = function() {
    View2D.dataInitialized = true;
    data.range = { max : -Infinity, min : Infinity};
    data.labelRange = {};
    data.values.forEach(function(d) {
      if (d > data.range.max) data.range.max = d;
      if (d < data.range.min) data.range.min = d;
      if (d > data.labelRange.max) data.labelRange.max = d;
      if (d < data.labelRange.min) data.labelRange.min = d;
    });
    View2D.drawData();
  }

  View2D.drawData = function() {
    /* value color legend */
    var formatToOne = d3.format(".1e");
    $('.min-value-text').text(formatToOne(data.range.min));
    $('.max-value-text').text(formatToOne(data.range.max));

    /* value color scale */
    var cs = d3.scaleLinear().domain([data.range.min, 0, data.range.max])
        .range([d3.rgb(PresetColor.blue), d3.rgb(PresetColor.white), d3.rgb(PresetColor.red)]);
    var valueColorScale = function(v) {
      return new THREE.Color(cs(v).toString());
    };

    // /* label color scale */
    // var cs2 = d3.scaleLinear().domain([labelRange.min, 0, labelRange.max])
    //     .range([d3.rgb(PresetColor.green), d3.rgb(PresetColor.black), d3.rgb(PresetColor.red)]);
    // var labelColorScale = function(l) {
    //   return new THREE.Color(cs2(l).toString());
    // };
    View2D.drawDataColor(valueColorScale);
  }

  View2D.drawDataColor = function(valueColorScale) {
    for (var i = 0; i < View2D.valueObject.geometry.faces.length; i ++) {
      var face = View2D.valueObject.geometry.faces[i];
      
      // update label color
      var l1 = data.labels[face.a];
      var l2 = data.labels[face.b];
      var l3 = data.labels[face.c];
      var labelColor = new THREE.Color(PresetColor.blackHex);
      if (l1 === l2 && l2 === l3) {
        // labelColor = labelColorScale(l1);
        if (l1 === 0) {
          labelColor = new THREE.Color(PresetColor.blackHex);
        }
        else if (l1 > 0) {
          labelColor = new THREE.Color(PresetColor.redHex);
        }
        else {
          labelColor = new THREE.Color(PresetColor.blueHex);
        }
      }
      View2D.labelObject.geometry.faces[i].color.set(labelColor);

      // update value color
      var v1 = data.values[face.a];
      var v2 = data.values[face.b];
      var v3 = data.values[face.c];
      View2D.valueObject.geometry.faces[i].vertexColors[0].set(valueColorScale(v1));
      View2D.valueObject.geometry.faces[i].vertexColors[1].set(valueColorScale(v2));
      View2D.valueObject.geometry.faces[i].vertexColors[2].set(valueColorScale(v3));
    }
    View2D.valueObject.geometry.colorsNeedUpdate = true;
    View2D.valueObject.material.needsUpdate = true;
    View2D.labelObject.geometry.colorsNeedUpdate = true;
    View2D.labelObject.material.needsUpdate = true;
    $('#loading').hide();
    $('.timestep-value').html(data.timestep);
    doneRendering = true;
  }

  View2D.updateRenderMethod = function(method) {
    if (method === 'label') {
      View2D.scene.remove(View2D.valueObject);
      View2D.scene.add(View2D.labelObject);
      $('.value-legend').hide();
      $('.label-legend').show();
      View2D.targetList = [View2D.labelObject];
    }
    else if (method === 'value'){
      View2D.scene.remove(View2D.labelObject);
      View2D.scene.add(View2D.valueObject);
      $('.value-legend').show();
      $('.label-legend').hide();
      View2D.targetList = [View2D.valueObject];
    }
  }

  View2D.updateRenderWireframe = function(renderWireframe) {
    View2D.valueObject.material.wireframe = renderWireframe;
    View2D.valueObject.geometry.colorsNeedUpdate = true;
    View2D.valueObject.material.needsUpdate = true;
  }

  function onDocumentMouseMove(event) {
    var isDialogShown = $('#connectDialog').hasClass('show');
    if (isDialogShown) return;
    if (data.values == undefined) return;
    var mouse = {};
    var parentBound = View2D.getView().getBoundingClientRect();
    var mouseX = event.clientX - parentBound.left;
    var mouseY = event.clientY - parentBound.top;
    // mouse.x = ( event.clientX / View2D.renderer.domElement.clientWidth ) * 2 - 1;
    // mouse.y = - ( event.clientY / View2D.renderer.domElement.clientHeight ) * 2 + 1;
    mouse.x = (mouseX / View2D.renderer.domElement.clientWidth) * 2 - 1;
    mouse.y = - (mouseY / View2D.renderer.domElement.clientHeight) * 2 + 1;
    View2D.raycaster.setFromCamera(mouse, View2D.camera);

    var intersects = View2D.raycaster.intersectObjects(View2D.targetList);

    if ( intersects.length > 0 ) {
      var face = intersects[0].face;  
      var show = hitFace(face, event);
      if (show) {
        $('.info-tooltip').show();
        $('.info-tooltip').css('left', event.clientX);
        $('.info-tooltip').css('top', event.clientY - parseFloat($('.info-tooltip').css('height')));
        $('#webGL-canvas').addClass('pointer-class');
        View2D.scene.add(View2D.circleObject);
      }
      else {
        $('.info-tooltip').hide();
        $('#webGL-canvas').removeClass('pointer-class');
        View2D.scene.remove(View2D.circleObject);
      }
    }
    else {
      $('.info-tooltip').hide();
      $('#webGL-canvas').removeClass('pointer-class');
      View2D.scene.remove(View2D.circleObject);
    }
  }

  function hitFace(face, event) {
    // calculate closest point 
    var parentBound = View2D.getView().getBoundingClientRect();
    var mouseX = event.clientX - parentBound.left;
    var mouseY = event.clientY - parentBound.top;
    var i1 = face.a;
    var i2 = face.b;
    var i3 = face.c;
    var c1 = [data.mesh.coords[i1*2] - 1.7, data.mesh.coords[i1*2+1]];
    var c2 = [data.mesh.coords[i2*2] - 1.7, data.mesh.coords[i2*2+1]];
    var c3 = [data.mesh.coords[i3*2] - 1.7, data.mesh.coords[i3*2+1]];
    var v1 = data.values[i1];
    var v2 = data.values[i2];
    var v3 = data.values[i3];
    var l1 = data.labels[i1];
    var l2 = data.labels[i2];
    var l3 = data.labels[i3];
    var squareC1 = (mouseX - c1[0])*(mouseX - c1[0]) + (mouseY - c1[1])*(mouseY - c1[1]);
    var squareC2 = (mouseX - c2[0])*(mouseX - c2[0]) + (mouseY - c2[1])*(mouseY - c2[1]);
    var squareC3 = (mouseX - c3[0])*(mouseX - c3[0]) + (mouseY - c3[1])*(mouseY - c3[1]);
    var i, c, v, l; 
    if (squareC1 > squareC2 && squareC1 > squareC3) {
      i = i1;
      c = c1;
      v = v1;
      l = l1;
    }
    else if (squareC2 > squareC1 && squareC2 > squareC3) {
      i = i2;
      c = c2;
      v = v2;
      l = l2;
    }
    else {
      i = i3;
      c = c3;
      v = v3;
      l = l3;
    }
    if (l == 0) {
      return false;
    }
    var formatToOne = d3.format(".1e");
    $('.info-indices')[0].innerHTML = i;
    $('.info-coords')[0].innerHTML = '(' + formatToOne(c[0]) + ', ' + formatToOne(c[1]) + ')';
    $('.info-values')[0].innerHTML = formatToOne(v);
    $('.info-labels')[0].innerHTML = l;
    View2D.circleObject.geometry = new THREE.CircleGeometry(.001, 32);
    View2D.circleObject.geometry.translate(c[0], c[1], 0);
    return true;
  }

  return View2D;
})();