var menuText = function() {
  this.dataName = "";
  this.frame = 0;
  this.resetTrackball = function () {
    cameraControls.reset();
  };
  // this.displayStats = false;
  this.renderMethod = 'value';
  this.renderWireframe = false;
  this.autoRefreshing = false;
  this.refreshInterval = 20;
  this.showTimestep2d = false;
  this.showLegend2d = false;
  this.showTimestep3d = false;
  this.showLegend3d = false;
  this.useSameTFEditor = false;
  this.enableAngle = false;
  this.startAngle = 0;
  this.endAngle = Math.PI * 2;

  this.reconnect = function () {
    connectionDialog.reconnect();
  };

  this.refresh = function() {
    requestData();
  };

  this.transferFunction = function() {
    ViewTF.toggle();
    $('.dg.main').parent().css('overflow-y', 'scroll');
  };
};

function initializeControlPanel(domElem) {
  var text = new menuText();
  var gui = new dat.GUI({autoPlace: false, width: 400});

  var f1 = gui.addFolder("2D Rendering");
  f1.add(text, 'renderMethod', ['value', 'label']).onChange(function() {
    View2D.updateRenderMethod(text.renderMethod);
  });
  f1.add(text, 'renderWireframe').onChange(function() {
    View2D.updateRenderWireframe(text.renderWireframe);
  });
  f1.add(text, 'showTimestep2d').onChange(function() {
    View2D.updateTimestepDisplay(text.showTimestep2d);
  });
  f1.add(text, 'showLegend2d').onChange(function() {
    View2D.updateLegendDisplay(text.showLegend2d);
  });
  f1.add(text, 'useSameTFEditor').onChange(function() {
    View2D.updateEnableSameTFEditor(text.useSameTFEditor);
  });
  f1.add(text, 'transferFunction');
  f1.open();

  var f2 = gui.addFolder("3D Rendering");
  // f2.add(text, "displayStats").onChange(function(val) {
  //   if (val) $("#stats").css({visibility: "visible"});
  //   else $("#stats").css({visibility: "hidden"});
  // });
  // f2.add(text, "saveImage");
  f2.add(text, 'showTimestep3d').onChange(function() {
    View3D.updateTimestepDisplay(text.showTimestep3d);
  });
  f2.add(text, 'showLegend3d').onChange(function() {
    View3D.updateLegendDisplay(text.showLegend3d);
  });
  f2.add(text, 'enableAngle').onChange(function() {
    View3D.updateEnableAngle();
  });
  f2.add(text, 'startAngle', 0, Math.PI * 2).onChange(function() {
    data.startAngle = text.startAngle;
    // TODO 
  });
  f2.add(text, 'endAngle', 0, Math.PI * 2).onChange(function() {
    data.endAngle = text.endAngle;
    // TODO
  });
  f2.open();

  var f3 = gui.addFolder("Connection");
  f3.add(text, 'refresh');
  f3.add(text, 'autoRefreshing').onChange(function() {
    if (text.autoRefreshing) {
      updateRepeatRequestingSpeed(text.refreshInterval);
    }
    else {
      cancelRepeatRequesting();
    }
  });
  f3.add(text, 'refreshInterval', 10, 300).onFinishChange(function() {
    if (text.autoRefreshing) {
      updateRepeatRequestingSpeed(text.refreshInterval);
    }
  });
  f3.add(text, 'reconnect');
  f3.open();
  domElem.appendChild(gui.domElement);
  return gui.domElement;
};

(function initialConnectDialog() {
  var clickConnect = function() {
    var ip = $('#hostAdress').val();
    var port = $('#port').val();
    connectToServer(ip, port);
    $('#connectDialog').modal('hide');
    $('#loading').show();
  };
  $('#connect').click(function() {
    clickConnect();
  });

  connectionDialog = {};

  connectionDialog.error = function() {
    $('#loading').hide();
    $('#connect-dialog-title')[0].innerHTML = 'Connection Error. Try Again';
    $('#close-dialog').hide();
    $('#close-dialog-cross').hide();
    $('#connectDialog').modal('show');
  };
  connectionDialog.closed = function() {
    $('#loading').hide();
    $('#connect-dialog-title')[0].innerHTML = 'Connection Closed. Try Again';
    $('#close-dialog').hide();
    $('#close-dialog-cross').hide();
    $('#connectDialog').modal('show');
  }
  connectionDialog.reconnect = function() {
    $('#loading').hide();
    $('#connect-dialog-title')[0].innerHTML = 'Reconnect to Server';
    $('#close-dialog').show();
    $('#close-dialog-cross').show();
    $('#connectDialog').modal('show');
  }

  document.querySelector('#connectDialog').addEventListener('keypress', function (e) {
      var key = e.which || e.keyCode;
      if (key === 13) { // 13 is enter
        console.log('connect');
        clickConnect();
      }
  });

  if (!DEBUG_MODE) $('#connectDialog').modal('show');
})();

// var stats = new Stats();
// stats.showPanel(0);
// stats.dom.id="stats";
// document.body.appendChild(stats.dom);
// $("#stats").css({visibility: "hidden"});

var layout = (function initialLayout() {
  var layout = {};
  var totalWidth = $(window).width();
  var left2Pct = (totalWidth - 400) / totalWidth / 2 * 100;
  var rightPct = 400 / totalWidth * 100;
  var config = {
      settings: {
        showPopoutIcon: false,
        showMaximiseIcon: false,
        showCloseIcon: false
      },
      content: [{
          type: 'row',
          content:[{
              type: 'column',
              content:[{
                  type: 'component',
                  componentName: '2D View',
                  isClosable: false,
                  componentState: {
                    cameraMenu: true,
                    viewId: '2d'
                  }
              },{
                  type: 'component',
                  componentName: 'FFT View',
                  isClosable: false,
                  componentState: {
                    // cameraMenu: true
                  }
              }],
              width: left2Pct
          },
          {
              type: 'column',
              content:[{
                  type: 'component',
                  componentName: '3D View',
                  isClosable: false,
                  componentState: {
                    cameraMenu: true,
                    viewId: '3d'
                  }
              },{
                  type: 'component',
                  componentName: 'Tree View',
                  isClosable: false,
                  componentState: {
                    cameraMenu: true,
                    viewId: 'tree'
                  }
              }],
              width: left2Pct
          },
          {
              type: 'column',
              content:[
              {
                  type: 'component',
                  componentName: 'Console',
                  componentState: {
                    cameraMenu: false
                  },
                  isClosable: false
              }],
              width: rightPct
          }]
      }]
  };
  
  var myLayout = new GoldenLayout(config);

  myLayout.registerComponent('2D View', function(container, componentState){
    layout.View2D = container.getElement()[0];
    View2D.initial();
    
    // var legendElem = $('.legend-div').detach();
    var legendElem = $('.legend-div').clone();
    legendElem.removeClass('template');
    legendElem.addClass('legend-div-2d');
    legendElem.find('linearGradient').attr('id', 'mainGradient2');
    $(layout.View2D).append(legendElem);
    
    // var timestepElem = $('.timestep-div').detach();
    var timestepElem = $('.timestep-div').clone();
    timestepElem.removeClass('template');
    timestepElem.addClass('timestep-div-2d');
    $(layout.View2D).append(timestepElem);
    
    new ResizeSensor($(layout.View2D), function(){ 
      View2D.resize();
    });
  });

  myLayout.registerComponent('3D View', function(container, componentState){
    layout.View3D = container.getElement()[0];
    View3D.initial();
    
    var legendElem = $('.legend-div').clone();
    legendElem.removeClass('template');
    legendElem.addClass('legend-div-3d');
    legendElem.find('linearGradient').attr('id', 'mainGradient3');
    $(layout.View3D).append(legendElem);

    var timestepElem = $('.timestep-div').clone();
    timestepElem.removeClass('template');
    timestepElem.addClass('timestep-div-3d');
    $(layout.View3D).append(timestepElem);
    
    new ResizeSensor($(layout.View3D), function(){ 
      View3D.resize();
    });
  });

  myLayout.registerComponent('Tree View', function(container, componentState){
    layout.ViewTree = container.getElement()[0];
    ViewTree.initial();
    requestTree();
    new ResizeSensor($(layout.ViewTree), function(){ 
      ViewTree.resize();
    });
  });

  myLayout.registerComponent('FFT View', function(container, componentState){
  });

  myLayout.registerComponent('Console', function(container, componentState){
    var guiDomElem = initializeControlPanel(container.getElement()[0]);
    window.guiDomElem = guiDomElem;
    $(guiDomElem).find('.close-button').remove();
    ViewTF.initial();
    var tfElem = $('#tf-holder').detach();
    $(guiDomElem).find('li.folder').first().append(tfElem);
  });

  myLayout.on('stackCreated', function(stack) {
    if (!stack.contentItems[0].config.componentState.cameraMenu) {
      return;
    }

    var dropdown = $($('template').html()),
      dropdownBtn = dropdown.find('.camera-menu-list');
    dropdownBtn.attr('viewId', stack.contentItems[0].config.componentState.viewId);

    var cameraAction = function(action, viewId) {
        var target;
        switch (viewId) {
          case '2d': 
            target = View2D;
            break;
          case '3d': 
            target = View3D;
            break;
          case 'tree': 
            target = ViewTree;
            break;
        }
        globalStatus.cameraActionTarget = target;
        switch (action) {
          case'save-image': 
            if (!data.mockLink) {
              data.mockLink = document.createElement("a");
            }
            if (viewId != '3d') { 
              data.mockLink.href = target.getRenderer().domElement.toDataURL('image/jpg');
            }
            else {
              data.mockLink.href = $('img.volren-image').attr('src');
            }
            data.mockLink.download = 'screenshot.png';
            data.mockLink.click();
            break;
          case 'save-camera':
            var cameraState = JSON.stringify(target.getCameraMatrix());
            var pom = document.createElement('a');
            pom.setAttribute('href', 'data:text/plain;charset=utf-8,' + encodeURIComponent(cameraState));
            pom.setAttribute('download', 'camera');
            if (document.createEvent) {
                var event = document.createEvent('MouseEvents');
                event.initEvent('click', true, true);
                pom.dispatchEvent(event);
            }
            else {
                pom.click();
            }
            break;
          case 'load-camera':
            $('.load-camera-file').click();
            break;
          case 'reset-camera':
            target.resetCamera();
            break;
        }
    };

    stack.header.controlsContainer.prepend(dropdown);

    dropdown.find('li').click(function(){
        cameraAction($(this).attr('action'), $(this).attr('viewid'))
    });

  });

  function readSingleFile(e) {
    var file = e.target.files[0];
    if (!file) {
      return;
    }
    var reader = new FileReader();
    reader.onload = function(e) {
      var contents = e.target.result;
      globalStatus.cameraActionTarget.setCameraMatrix(contents);
    };
    reader.readAsText(file);
  }
  document.getElementsByClassName('load-camera-file')[0]
    .addEventListener('click', function() {
      this.value = null;
    });
  document.getElementsByClassName('load-camera-file')[0]
    .addEventListener('change', readSingleFile, false);

  myLayout.init();
  return layout;
})();

