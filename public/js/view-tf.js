var ViewTF = (function() {
  var ViewTF = {};

  ViewTF.getView = function() {
    return layout.ViewTF;
  };
  
  ViewTF.initial = function() {
    $('#tf-holder').tfWidget(callback);
    function callback(controlPoints, tfArray) {
      requestImage(tfArray);
    }
  };

  return ViewTF;
})();