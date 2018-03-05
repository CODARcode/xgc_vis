var ViewControl = (function() {
  var ViewControl = {};

  ViewControl.getView = function() {
    return layout.ViewControl;
  };

  ViewControl.initial = function() {
    var guiDomElem = initializeControlPanel(ViewControl.getView());
    $(guiDomElem).find('.close-button').remove();
    ViewTF.initial();
    var tfElem = $('#tf-holder').detach();
    $(guiDomElem).find('li.folder').first().append(tfElem);
    $('.dg.main').parent().css('overflow-y', 'scroll');
  };

  return ViewControl;
})();
