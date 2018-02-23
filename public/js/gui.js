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

  this.saveImage = function() {
    window.open( renderer.domElement.toDataURL( 'image/png' ), 'screenshot' );
  };

  this.reconnect = function () {
    connectionDialog.reconnect();
  };

  this.refresh = function() {
    requestData();
  };

  this.resetCamera = function() {
    controls.reset();
  };

  this.transferFunction = function() {
    ViewTF.toggle();
  }
};

function initializeControlPanel(domElem) {
  var text = new menuText();
  var gui = new dat.GUI({autoPlace: false, width: 400});

  // var f1 = gui.addFolder("Simulation");
  // f1.add(text, 'dataName');
  // f1.add(text, 'timestep'); 
  // f1.open();

  var f2 = gui.addFolder("3D Rendering");
  // f2.add(text, "displayStats").onChange(function(val) {
  //   if (val) $("#stats").css({visibility: "visible"});
  //   else $("#stats").css({visibility: "hidden"});
  // });
  // f2.add(text, "saveImage");
  f2.add(text, 'renderMethod', ['value', 'label']).onChange(function() {
    View2D.updateRenderMethod(text.renderMethod);
  });
  f2.add(text, 'renderWireframe').onChange(function() {
    View2D.updateRenderWireframe(text.renderWireframe);
  });
  f2.add(text, 'transferFunction');
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
                  componentState: {}
              },{
                  type: 'component',
                  componentName: 'FFT View',
                  isClosable: false,
                  componentState: {}
              }],
              width: left2Pct
          },
          {
              type: 'column',
              content:[{
                  type: 'component',
                  componentName: '3D View',
                  isClosable: false,
                  componentState: {}
              },{
                  type: 'component',
                  componentName: 'Tree View',
                  isClosable: false,
                  componentState: {}
              }],
              width: left2Pct
          },
          {
              type: 'column',
              content:[
              {
                  type: 'component',
                  componentName: 'Console',
                  componentState: {},
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
    new ResizeSensor($(layout.View2D), function(){ 
      View2D.resize();
    });
  });

  myLayout.registerComponent('3D View', function(container, componentState){
    layout.View3D = container.getElement()[0];
    View3D.initial();
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
    // var elem = container.getElement()[0];
    // var legendElem = $('#legend-div').detach();
    // $(elem).append(legendElem);
    var tfElem = $('#tf-holder').detach();
    // $(elem).append(tfElem);
    $(guiDomElem).find('li.folder').first().append(tfElem);
    // layout.ViewTF = elem;
  });

  // myLayout.on( 'stackCreated', function( stack ){

  //   //HTML for the colorDropdown is stored in a template tag
  //   var colorDropdown = $( $( 'template' ).html() ),
  //       colorDropdownBtn = colorDropdown.find( '.selectedColor' );


  //   var setColor = function( color ){
  //       var container = stack.getActiveContentItem().container;

  //       // Set the color on both the dropDown and the background
  //       colorDropdownBtn.css( 'background-color', color );
  //       container.getElement().css( 'background-color', color );

  //       // Update the state
  //       container.extendState({ color: color });
  //   };

  //   // Add the colorDropdown to the header
  //   stack.header.controlsContainer.prepend( colorDropdown );

  //   // Update the color initially and whenever the tab changes
  //   stack.on( 'activeContentItemChanged', function( contentItem ){
  //        setColor( contentItem.container.getState().color );
  //   });
       
  //   // Update the color when the user selects a different color
  //   // from the dropdown
  //   colorDropdown.find( 'li' ).click(function(){
  //       setColor( $(this).css( 'background-color' ) );
  //   });
  // });

  myLayout.init();
  return layout;
})();