var ViewAngle = (function() {
  var ViewAngle = {};

  ViewAngle.getView = function() {
    return layout.ViewAngle;
  };

  ViewAngle.draw = function() {
    ViewAngle.initialized = true;
    var elem = ViewAngle.getView();
    var width = $(elem).width();
    var height = $(elem).height();
    var abs = Math.max(Math.abs(data.angleValueMin), Math.abs(data.angleValueMax));
    ViewAngle.chart = Highcharts.chart('angle-chart-container', {
      title:{
        text:''
      },
      legend: {
          enabled: false
      },
      series: [
        { 
          name: 'value',
          data: data.angleData
        }
      ],
      chart: {
        type: 'area'
      },
      yAxis: {
        min: -abs,
        max: abs
      }
    });
    ViewAngle.chart.setSize(width, height);
  };

  ViewAngle._draw = function() {
    
    var elem = ViewAngle.getView();
    var width = $(elem).width();
    var height = $(elem).height();
    $('.angle-svg').attr('width', width);
    $('.angle-svg').attr('height', height);
    if (!ViewAngle.initialized) {
      var path = $('<g><path class="line"></path></g>');
      var xAxis = $('<g class="line-x-axis"></g>');
      var yAxis = $('<g class="line-y-axis"></g>');
      $('.angle-svg').append(path);
      $('.angle-svg').append(xAxis);
      $('.angle-svg').append(yAxis);
      ViewAngle.initialized = true;
    }
    ViewAngle.drawLine();

    
  }

  ViewAngle.drawLine = function() {
    var elem = ViewAngle.getView();
    var width = $(elem).width();
    var height = $(elem).height();
    var xGap = 30;
    var yGap = 20;
    var xStart = xGap;
    var xEnd = width - xGap;
    var yStart = yGap;
    var yEnd = height - yGap;
    var xScale = d3.scaleLinear().domain([data.thetaMin, data.thetaMax])
        .range([xStart, xEnd]);
    var yScale = d3.scaleLinear().domain([data.angleValueMin, data.angleValueMax])
        .range([yEnd, yStart]);
    var area = d3.area()
      .x(function(d) {
        var x = xScale(d.theta);
        if (x == 0) {
          console.log(d);
        }
        return x;
      })
      .y0(yScale(0))
      .y1(function(d) {
        var y = yScale(d.value); 
        if (y == 0) {
          console.log(d);
        }
        return y;
      });
    d3.select('svg.angle-svg')
      .select('.line')
      .datum(data.angleData)
      .attr("fill", PresetColor.blue)
      .attr('stroke', 'none')
      .attr("d", area);
    
    // x axis
    d3.select('g.line-x-axis')
      .attr("transform", "translate(0," + yScale(0) + ")")
      .call(d3.axisBottom(xScale));

    // y axis
    d3.select('g.line-y-axis')
      .attr("transform", "translate(" + xStart + ",0)")
      .call(d3.axisLeft(yScale));

    $('.angle-svg').html($('.angle-svg').html());


  };

  ViewAngle.resize = function() {
    if (!ViewAngle.initialized) return;
    var elem = ViewAngle.getView();
    var width = $(elem).width();
    var height = $(elem).height();
    ViewAngle.chart.setSize(width, height);
    ViewAngle.chart.reflow();
  };

  return ViewAngle;
})();