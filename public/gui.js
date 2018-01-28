var menuText = function() {
  this.dataName = "";
  this.frame = 0;
  this.resetTrackball = function () {
    cameraControls.reset();
  };
  this.displayStats = false;
  this.renderMethod = 'label';
  this.renderWireframe = false;
  this.autoRefreshing = false;
  this.refreshInterval = 20;

  this.saveImage = function() {
    window.open( renderer.domElement.toDataURL( 'image/png' ), 'screenshot' );
  };

  this.saveTrackball = function () {
    var trac = {
      target: cameraControls.target,
      position: cameraControls.object.position,
      up: cameraControls.object.up
    };
    var blob = new Blob([JSON.stringify(trac)], {type: "text/plain;charset=utf-8"});
    saveAs(blob, "camera.json");
  };

  this.loadTrackball = function () {
    if (!window.File || !window.FileReader || !window.FileList || !window.Blob) {
      alert('The File APIs are not fully supported in this browser.');
      return;
    }

    var evt = document.createEvent("MouseEvents");
    evt.initEvent("click", true, false);
    file_open.dispatchEvent(evt);
    file_open.onchange = function() {
      var path = file_open.value;
      var f = file_open.files[0];
      var reader = new FileReader();
      
      reader.onload = function(evt) {
        var str = evt.target.result;
        var trac = JSON.parse(str);
        cameraControls.load(trac.target, trac.position, trac.up);
      }

      reader.readAsText(f);
    }
  };

  this.reconnect = function () {
    connectionDialog.reconnect();
  };
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
  f2.add(text, 'renderMethod', ['label', 'value']).onChange(function() {
    updateRenderMethod(text.renderMethod);
  });
  f2.add(text, 'renderWireframe').onChange(function() {
    updateRenderWireframe(text.renderWireframe);
  });
  f2.open();

  var f3 = gui.addFolder("Connection");
  f3.add(text, 'autoRefreshing').onChange(function() {
    if (text.autoRefreshing) {
      updateRepeatRequestingSpeed(text.repeatInterval);
    }
    else {
      cancelRepeatRequesting();
    }
  });
  f3.add(text, 'refreshInterval', 10, 300).onFinishChange(function() {
    updateRepeatRequestingSpeed(text.repeatInterval);
  });
  f3.add(text, 'reconnect');
  f3.open();

  var f4 = gui.addFolder("Trackball");
  f4.add(text, "saveTrackball");
  f4.add(text, "loadTrackball");
  f4.add(text, "resetTrackball");
};

(function initialConnectDialog() {
  $('#connect').click(function() {
    var ip = $('#hostAdress').val();
    var port = $('#port').val();
    connectToServer(ip, port);
    $('#connectDialog').modal('hide');
    $('#loading').show();
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

  $('#connectDialog').modal('show');
})();