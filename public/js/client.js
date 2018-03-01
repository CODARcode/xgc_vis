var ws;
var Client = {};

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
    // console.log(returnTreeData);
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
    data.treeDepth = 0;
    var sizeScale = 500;
    while (stack.length > 0) {
      var curNode = stack.pop();
      data.treeNodeCount ++;
      if (data.treeDepth < curNode.depth) data.treeDepth = curNode.depth;
      var x = curNode.x;
      var y = curNode.y;
      curNode.x = y * Math.cos(x) * sizeScale;
      curNode.y = y * Math.sin(x) * sizeScale;
      data.treeData[curNode.id] = curNode;

      if (curNode.data.extremum_val > extremumMax) extremumMax = curNode.data.extremum_val;
      if (curNode.data.extremum_val < extremumMin) extremumMin = curNode.data.extremum_val;
      if (curNode.data.saddle_val > saddleMax) saddleMax = curNode.data.saddle_val;
      if (curNode.data.saddle_val < saddleMin) saddleMin = curNode.data.saddle_val;
      if (curNode.children) {
        stack = stack.concat(curNode.children);
      }
    }
    console.log(extremumMin, extremumMax);
    console.log(saddleMin, saddleMax);
    data.treeMinMax = {
      extremumMax: extremumMax,
      extremumMin: extremumMin,
      saddleMin: saddleMin,
      saddleMax: saddleMax
    };
    t2 = new Date();
    // console.log(t2.getTime() - t1.getTime());

    var maxExtremumAbs = Math.max(Math.abs(data.treeMinMax.extremumMin), data.treeMinMax.extremumMax);
    var maxSaddleAbs = Math.max(Math.abs(data.treeMinMax.saddleMin), data.treeMinMax.saddleMax);
    var maxAbs = Math.max(maxExtremumAbs, maxSaddleAbs);
    var maxNodeHeight = 1 * sizeScale;
    var heightScale = d3.scaleLinear().domain([0, maxAbs])
        .range([0, maxNodeHeight]);
    for (var id in data.treeData) {
      var extremumY = heightScale(data.treeData[id].data.extremum_val);
      var saddleY = heightScale(data.treeData[id].data.saddle_val);
      data.treeData[id].saddleY = saddleY;
      data.treeData[id].extremumY = extremumY;
      if (isNaN(saddleY) || isNaN(extremumY)) {
        debugger;
        console.log(data.treeData[id]);
      }
    }

    ViewTree.drawRadialTree();
  });
}

function requestData() {
  if (VIEW2D_OFF) return;
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
  if (VIEW3D_OFF) return;
  console.log("requesting multiple slice raw data...");
  var msg = {
    type: "requestMultipleSliceRawData"
  };

  if (ws.readyState == 1) ws.send(JSON.stringify(msg));
  else connectToServer();
}

function requestImageWait(tfArray) {
  if (volrenTimer) {
    clearTimeout(volrenTimer);
  }
  volrenTimer = setTimeout(function() {
    requestImage(tfArray);
  }, 1000 * VOLREN_TIME);
}

function requestImage(tfArray) {
  if (IMAGE_OFF) return;
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
  if (!tfArray) {
    tfArray = ViewTF.getTF();
  }
  if (tfArray && tfArray.length === 1024) {
    var containsNaN = false;
    for (var i = 0; i < tfArray.length; i ++) {
      if (isNaN(tfArray[i])) {
        containsNaN = true;
        break;
      }
    }
    if (!containsNaN) {
      msg.tf = tfArray;
    }
  }
  console.log('POST! ', msg);
  if (DEBUG_MODE || VIEW3D_OFF) return;

  if (ws && ws.readyState == 1) {
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
  }

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
  requestImage();
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
  // console.log('is binary: ', isBinary);
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
      // requestImage();
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

