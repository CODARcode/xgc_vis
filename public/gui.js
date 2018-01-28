var menuText = function() {
  this.dataName = "";
  this.frame = 0;
  this.resetTrackball = function () {
    cameraControls.reset();
  };
  this.displayStats = false;

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
  f2.add(text, "saveImage");
  f2.open();

  var f3 = gui.addFolder("Trackball");
  f3.add(text, "saveTrackball");
  f3.add(text, "loadTrackball");
  f3.add(text, "resetTrackball");
};
