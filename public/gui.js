var menuText = function() {
  this.dataName = "";
  this.frame = 0;
  this.resetTrackball = function () {
    cameraControls.reset();
  };
  this.displayStats = false;
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
  }

  this.resetCamera = function() {
    controls.reset();
  }
};

function initializeControlPanel () {
  var text = new menuText();
  var gui = new dat.GUI();

  // var f1 = gui.addFolder("Simulation");
  // f1.add(text, 'dataName');
  // f1.add(text, 'timestep'); 
  // f1.open();

  var f2 = gui.addFolder("3D Rendering");
  f2.add(text, "displayStats").onChange(function(val) {
    if (val) $("#stats").css({visibility: "visible"});
    else $("#stats").css({visibility: "hidden"});
  });
  // f2.add(text, "saveImage");
  f2.add(text, 'renderMethod', ['value', 'label']).onChange(function() {
    updateRenderMethod(text.renderMethod);
  });
  f2.add(text, 'renderWireframe').onChange(function() {
    updateRenderWireframe(text.renderWireframe);
  });
  f2.add(text, 'resetCamera');
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

  $('#connectDialog').modal('show');
})();

