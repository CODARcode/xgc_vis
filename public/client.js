const wsUri = "ws://localhost:9002";
var ws;

function requestMesh() {
  console.log("requesting mesh");
  var msg = {
    type: "requestMesh"
  };

  if (ws.readyState == 1) ws.send(JSON.stringify(msg));
  else connectToServer();
}

function connectToServer() {
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
  console.log("connected to server");
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
  console.log(msg);

  /* 
  if (msg.type == "dbList") {
    updateDBList(msg.data);
  }
  else if (msg.type == "dataInfo") {
    updateDataInfo(msg.dataInfo, msg.events);
  }
  else if (msg.type == "vlines") {
    updateVlines(msg.data.vlines);
    updateDistances(msg.data.dist);
  } */
}

function onError(evt)
{
  console.log("error");
}

connectToServer();
