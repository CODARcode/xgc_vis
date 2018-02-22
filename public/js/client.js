var ws;

function requestMesh() {
  console.log("requesting mesh...");
  var msg = {
    type: "requestMesh"
  };

  if (ws.readyState == 1) ws.send(JSON.stringify(msg));
  else connectToServer();
}

function requestTree() {
  d3.json('data/branches.00001.json', function(returnTreeData) {
    console.log(returnTreeData);
    var stratify = d3.stratify()
      .parentId(function(d) { return d.parent; })
      .id(function(d) {return d.id;});
    var tree = d3.tree()
        .size([2 * Math.PI, 1])
        .separation(function(a, b) { return (a.parent == b.parent ? 1 : 2) / a.depth; });
    var t1 = new Date();
    data.d3treeRoot = tree(stratify(returnTreeData));
    var t2 = new Date();
    // console.log(t2.getTime() - t1.getTime());
    
    var extremumMin = saddleMin = +Infinity;
    var extremumMax = saddleMax = -Infinity;

    var stack = [data.d3treeRoot];
    t1 = new Date();
    data.treeNodeCount = 0;
    while (stack.length > 0) {
      var curNode = stack.pop();
      data.treeNodeCount ++;
      var x = curNode.x;
      var y = curNode.y;
      curNode.x = y * Math.cos(x);
      curNode.y = y * Math.sin(x);
      data.treeData[curNode.id] = curNode;

      if (curNode.extremum_val > extremumMax) extremumMax = curNode.extremum_val;
      if (curNode.extremum_val < extremumMin) extremumMin = curNode.extremum_val;
      if (curNode.saddle_val > saddleMax) saddleMax = curNode.saddle_val;
      if (curNode.saddle_val < saddleMin) saddleMin = curNode.saddle_val;
      if (curNode.children) {
        stack = stack.concat(curNode.children);
      }
    }
    t2 = new Date();
    // console.log(t2.getTime() - t1.getTime());
    ViewTree.drawRadialTree();
  });
}

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

function requestSingleSliceRawData() {
  console.log("requesting single slice raw data...");
  var msg = {
    type: "requestSingleSliceRawData"
  };

  if (ws.readyState == 1) ws.send(JSON.stringify(msg));
  else connectToServer();
}

function requestMultipleSliceRawData() {
  console.log("requesting multiple slice raw data...");
  var msg = {
    type: "requestMultipleSliceRawData"
  };

  if (ws.readyState == 1) ws.send(JSON.stringify(msg));
  else connectToServer();
}

function requestImage() {
  var wh = View3D.getWH();
  var msg = {
    type: "requestVolren", 
    matrix: View3D.getMatrix(),
    fov: View3D.getCameraFov(),
    width: wh.width,
    height: wh.height,
    near: View3D.getCameraNear(),
    far: View3D.getCameraFar()
  };
  console.log(msg);

  var url = 'http://' + ws.url.substr(5);
  fetch(url + 'requestVolren', {
    body: JSON.stringify(msg),
    method: 'POST'
  })
  .then(function(response) {
    return response.blob();
  }).then(function(imageDataBlob) {
    View3D.updateImage(imageDataBlob);
  });

  // if (ws.readyState == 1) ws.send(JSON.stringify(msg));
  // else connectToServer();
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
  
  ws.binaryType = "arraybuffer";
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
  requestTree();
}

function onClose(evt)
{
  console.log("connection closed");
  connectionDialog.error();
}

function onMessage(evt)
{
  console.log(evt);
  var isBinary = !(typeof evt.data === 'string');
  console.log('is binary: ', isBinary);
  if (!isBinary) {
    var msg = JSON.parse(evt.data);
    if (msg.type == 'mesh') {
      console.log(msg.data);
      data.mesh = msg.data;
      data.nPhi = msg.data.nPhi;
      data.nNodes = msg.data.nNodes;
      View2D.drawMesh();
      requestData();
      View3D.drawMesh();
      requestImage();
      requestMultipleSliceRawData();
    } else if (msg.type == 'data') {
      console.log(msg);
      data.values = msg.data;
      data.labels = msg.labels;
      data.timestep = msg.timestep;
      View2D.updateData();
    }
  }
  else {
    var dataView = new DataView(evt.data);
    var dataType = dataView.getInt32(0, true);
    console.log(dataType);
    console.log(evt.data);
    if (dataType === 10) { // single slice raw data
      var singleSliceArray = new Float32Array(evt.data, 4);
      console.log(singleSliceArray);
      data.singleSliceArray = singleSliceArray;
    }
    else if (dataType === 11) { // multiple slice raw data
      data.slices = [];
      for (var i = 0; i < data.nPhi; i ++) {
        var sliceArray = new Float32Array(evt.data, 4 * (1 + data.nNodes * i), data.nNodes);
        data.slices.push(sliceArray);
      }
      console.log(data.slices);
      View3D.updateData();
    }
    else if (dataType === 13) {
      var imageArray = new Uint8Array(evt.data, 4);
      View3D.updateImage(imageArray);
    }
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

