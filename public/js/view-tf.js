var ViewTF = (function() {
  var ViewTF = {};

  ViewTF.getView = function() {
    return layout.ViewTF;
  };
  
  ViewTF.initial = function() {
    window.tfWidget = $('#tf-holder').tfWidget(callback);
    console.log(tfWidget);
    function callback(controlPoints, tfArray) {
      data.tfArray = tfArray;
      requestImageWait(tfArray);
    }
  };

  ViewTF.getTF = function() {
    return data.tfArray;
  };

  return ViewTF;
})();