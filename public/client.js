const wsUri = "ws://localhost:9002";
var ws;

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

function connectToServer() {
  console.log("connecting to server...");
  // ws = new WebSocket("ws://red.mcs.anl.gov:8080");
  ws = new WebSocket(wsUri);
  
  // ws.binaryType = "arraybuffer";
  ws.onopen = onOpen;
  ws.onclose = onClose;
  ws.onerror = onError;
  ws.onmessage = onMessage;
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
}

function onMessage(evt)
{
  var msg = JSON.parse(evt.data);
  if (msg.type == "mesh") {
    updateMesh(msg.data);
  } else if (msg.type == "data") {
    console.log(msg.data);
  }
}

function onError(evt)
{
  console.log("error");
}

connectToServer();
