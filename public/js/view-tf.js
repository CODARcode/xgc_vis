var ViewTF = (function() {
  var ViewTF = {};

  ViewTF.toggle = function() {
    $('#tf-holder').toggle();
  };
  
  ViewTF.initial = function() {
    ViewTF.initialized = true;
    var settings = {
      controlPoints: [
        {
          percent: 0, 
          alpha: 1, 
          rgb: PresetColor.blue
        },
        {
          percent: .5, 
          alpha: 0, 
          rgb: PresetColor.white
        },
        {
          percent: 1,
          alpha: 1,
          rgb: PresetColor.red
        }
      ]
    };
    var tfWidget = $('#tf-holder').tfWidget(callback, settings);
    function callback(controlPoints, tfArray) {
      data.tfArray = tfArray;
      requestImageWait(tfArray);
    }
    tfWidget.hide();
  };

  ViewTF.getTF = function() {
    return data.tfArray;
  };

  return ViewTF;
})();