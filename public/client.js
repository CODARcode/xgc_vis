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

function requestData() {
  console.log("requesting data...");
  var msg = {
    type: "requestData"
  };
  
  if (ws.readyState == 1) ws.send(JSON.stringify(msg));
  else connectToServer();
}

function connectToServer(ip, port) {
  console.log("connecting to server...");
  // ws = new WebSocket("ws://red.mcs.anl.gov:8080");
  uri = 'ws://' + ip + ':' + port;
  ws = new WebSocket(uri);
  
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
    console.log(msg.data);
    console.log(msg.labels);
    updateData(msg.data, msg.labels);
  }
}

function onError(evt)
{
  connectionDialog.error();
}

(function repeatRequestingData() {
  function repeatRequest() {
    requestData();
  }
  var repeatRequesting;

  updateRepeatRequestingSpeed = function(intervalSeconds) {
    repeatRequest();
    if (repeatRequesting != undefined) clearInterval(repeatRequesting);
    repeatRequesting = setInterval(repeatRequest, intervalSeconds*1000);
  };

  cancelRepeatRequesting = function() {
    if (repeatRequesting != undefined) clearInterval(repeatRequesting);
  };
})();

