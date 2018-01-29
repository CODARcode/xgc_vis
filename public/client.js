const wsUri = "ws://192.168.0.14:9003";
var ws;
var data = {};

function requestMesh() {
  console.log("requesting mesh...");
  var msg = {
    type: "requestMesh"
  };

  if (ws.readyState == 1) ws.send(JSON.stringify(msg));
  else connectToServer();
}

var doneRendering = true;
function requestData() {
  if (!doneRendering) return;
  doneRendering = false;
  console.log("requesting data...");
  var msg = {
    type: "requestData"
  };
  if (data.timestep != undefined) 
    msg['client_current_timestep'] = data.timestep;
  
  if (ws.readyState == 1) ws.send(JSON.stringify(msg));
  else connectToServer();
}

function connectToServer(ip, port) {
  console.log("connecting to server...");
  // ws = new WebSocket("ws://red.mcs.anl.gov:8080");
  uri = 'ws://' + ip + ':' + port;
  ws = new WebSocket(uri);

  setTimeout(function () {
    if (ws.readyState != 0 && ws.readyState != 1) {
        ws.onerror();
    }
  }, 3000);
  
  // ws.binaryType = "arraybuffer";
  ws.onopen = onOpen;
  ws.onclose = onClose;
  ws.onerror = onError;
  ws.onmessage = onMessage;
  $('#loading').show();
}

function onOpen(evt)
{
  console.log("connected to server.");
  requestMesh();
  // requestFrame(currentFrame);
}

function onClose(evt)
{
  console.log("connection closed");
  connectionDialog.error();
}

function onMessage(evt)
{
  var msg = JSON.parse(evt.data);
  if (msg.type == "mesh") {
    console.log(msg.data);
    updateMesh(msg.data);
  } else if (msg.type == "data") {
    console.log(msg);
    updateData(msg.data, msg.labels, msg.timestep);
  }
}

function onError(evt)
{
  connectionDialog.error();
}

function isDialogShown() {
  return $('#connectDialog').hasClass('show');
}

globalStatus = {};
(function repeatRequestingData() {
  function repeatRequest() {
    if (isDialogShown()) return;
    requestData();
  }
  var repeatRequesting;

  updateRepeatRequestingSpeed = function(intervalSeconds) {
    if (!Number.isInteger(intervalSeconds)) return;
    repeatRequest();
    if (repeatRequesting != undefined) clearInterval(repeatRequesting);
    repeatRequesting = setInterval(repeatRequest, intervalSeconds*1000);
  };

  cancelRepeatRequesting = function() {
    if (repeatRequesting != undefined) clearInterval(repeatRequesting);
  };
})();

