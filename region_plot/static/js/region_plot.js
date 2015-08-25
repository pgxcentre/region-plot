
//////////////////////////
// Interesting function //
//////////////////////////

var translate = function(x, y) {
  return "translate(" + x + ", " + y + ")";
}


////////////////////////////////////
// Parameters for the region plot //
////////////////////////////////////

var lowerOpacity = 0.08;

var regionChrom = {{ chrom }},
    regionStart = {{ start }} / 1e6,
    regionEnd = {{ end }} / 1e6;

var margin = { left: 45, top: 20, right: 60, bottom: 45 },
    width = $(window).width() - margin.left - margin.right,
    height = width * 1/3;

var pointSize = 3,
    diamondSize = 30;


//////////////////////////
// The scatter plot SVG //
//////////////////////////

var svg = d3.select("#scatter").append("svg")
  .attr("width", (width + margin.left + margin.right))
  .attr("height", (height + margin.top + margin.bottom));

var wrapper = svg.append("g")
  .attr("class", "scatterWrapper")
  .attr("transform", translate(margin.left, margin.top))


/////////////////////////////////////
// The X axis and scale (position) //
/////////////////////////////////////

var xScale = d3.scale.linear()
  .range([0, width])
  .domain([regionStart, regionEnd])
  .nice();

var xAxis = d3.svg.axis()
  .orient("bottom")
  .scale(xScale);

wrapper.append("g")
  .attr("class", "x axis")
  .attr("transform", translate(5, height))
  .call(xAxis);

wrapper.append("g")
  .append("text")
  .attr("class", "x label")
  .attr("text-anchor", "middle")
  .attr("transform", translate(width / 2, height + margin.bottom - 10))
  .text("Position on chr" + regionChrom + " (Mb)");


//////////////////////////////////
// The Y axis and scale (assoc) //
//////////////////////////////////

var assocScale = d3.scale.linear()
  .range([height, 0])
  .domain([-0.2, {{ assoc_max }}]);

// The y axis (assoc)
var assocAxis = d3.svg.axis()
  .orient("left")
  .scale(assocScale);

// Append the y-axis (assoc
wrapper.append("g")
  .attr("class", "y axis")
  .attr("transform", translate(5, -5))
  .call(assocAxis);

wrapper.append("g")
  .append("text")
  .attr("class", "y label")
  .attr("text-anchor", "middle")
  .attr("transform", translate(-margin.left/2, height/2) + " rotate(-90)")
  .text("-log(p)");


/////////////////////
// The actual data //
/////////////////////

// Imputed data
var imputedGroup = wrapper.append("g")
  .attr("class", "imputedWrapper");

imputedGroup.selectAll("imputedPoint").data(imputedData)
  .enter().append("path")
  .attr("class", function(d) { return "point " + d.range + " imputed " + d.snp; })
  .attr("d", d3.svg.symbol().type("diamond").size(30))
  .attr("transform", function(d) {
    return translate(xScale(d.pos / 1e6), assocScale(d.assoc));
  })
  .on("mouseover", showTooltip)
  .on("mouseout", removeTooltip);

// Genotyped data
var genotypedGroup = wrapper.append("g")
  .attr("class", "genotypedWrapper");

genotypedGroup.selectAll("genotypedPoint").data(genotypedData)
  .enter().append("circle")
    .attr("class", function(d) { return "point " + d.range + " genotyped " + d.snp; })
    .attr("cx", function(d) { return xScale(d.pos / 1e6); })
    .attr("cy", function(d) { return assocScale(d.assoc); })
    .attr("r", 3)
    .on("mouseover", showTooltip)
    .on("mouseout", removeTooltip);


////////////////////////////////
// The tooltips for the point //
////////////////////////////////

function showTooltip(d) {

  var element = d3.selectAll(".point." + d.snp);

  $(element).popover({
    placement: "auto top",
    container: "#scatter",
    trigger: "manual",
    html: true,
    content: tooltipContent(d),
  });

  if (d.type === "genotyped") {
    element.transition().duration(200).attr("r", pointSize * 3);
  } else {
    element.transition().duration(200).attr("d", d3.svg.symbol().type("diamond").size(diamondSize * 12));
  }

  $(element).popover("show");

}

function removeTooltip(d) {

  var element = d3.selectAll(".point." + d.snp);

  if (d.type === "genotyped") {
    element.transition().duration(200).attr("r", pointSize);
  } else {
    element.transition().duration(200).attr("d", d3.svg.symbol().type("diamond").size(diamondSize));
  }

  $(".popover").each(function() { $(this).remove(); });

}

function tooltipContent(d) {

  var name = "<span style='font-size: 11px; text-align: center;'>" + d.snp + " (" + d.type + ")</span>";
  var assoc = "<span style='font-size: 11px; text-align: center;'>p = " + d.p + "</span>";
  var corr = "<span style='font-size: 11px; text-align: center;'>r2 = " + d.ld + "</span>";

  return name + "<br>" + assoc + "<br>" + corr;

}


////////////////
// The legend //
////////////////

var legendMargin = { left: 5, top: 10, right: 5, bottom: 10 },
    legendWidth = 145,
    legendHeight = 145;

var svgLegend = d3.select("#scatterlegend")
  .style("top", margin.top + "px")
  .style("right", margin.right + "px")
  .style("width", (legendWidth + legendMargin.left + legendMargin.right) + "px")
  .style("height", (legendHeight + legendMargin.top + legendMargin.bottom) + "px")
  .append("svg")
    .attr("width", (legendWidth + legendMargin.left + legendMargin.right))
    .attr("height", (legendHeight + legendMargin.top + legendMargin.bottom));

var legendWrapper = svgLegend.append("g")
  .attr("class", "legendWrapper")
  .attr("transform", translate(legendMargin.left, legendMargin.top));


////////////////////////
// The legend content //
////////////////////////

var rectSize = 15,
    rowHeight = 20,
    maxWidth = 144,
    rectValues = [
      { label: "r2 <= 0.2", range: "lower_02" },
      { label: "r2 <= 0.4", range: "lower_04" },
      { label: "r2 <= 0.6", range: "lower_06" },
      { label: "r2 <= 0.8", range: "lower_08" },
      { label: "r2 <= 1.0", range: "lower_1" },
    ]
    shapeValues = [
      { label: "Imputed", shape: "diamond", type: "imputed" },
      { label: "Genotyped", shape: "circle", type: "genotyped" },
    ];

var pointLegend = legendWrapper.selectAll("legendShape")
  .data(shapeValues).enter().append("g")
  .attr("class", "legendShape")
  .attr("transform", function(d, i) { return translate(rectSize / 2, i * rowHeight); })
  .style("cursor", "pointer")
  .on("mouseover", selectShapeLegend(lowerOpacity))
  .on("mouseout", selectShapeLegend(1))
  .on("click", clickShapeLegend);

pointLegend.append("path")
  .attr("class", function(d) { return "legend shape " + d.type; })
  .attr("d", d3.svg.symbol().type(function(d) { return d.shape; }).size(diamondSize * 3))

pointLegend.append("text")
  .attr("transform", translate(rectSize, 0))
  .attr("class", "legendText")
  .style("font-size", "10px")
  .style("fill", "#000000")
  .attr("dy", ".35em")
  .text(function(d) { return d.label; });

var legend = legendWrapper.selectAll(".legendSquare")
  .data(rectValues).enter().append("g")
    .attr("class", "legendSquare")
    .attr("transform", function(d, i) { return translate(0, (i+2) * rowHeight); })
    .style("cursor", "pointer")
    .on("mouseover", selectLegend(lowerOpacity))
    .on("mouseout", selectLegend(1))
    .on("click", clickLegend);

legend.append("rect")
  .attr("class", function(d) { return "legend rectangle " + d.range; })
  .attr("width", rectSize)
  .attr("height", rectSize)

legend.append("text")
  .attr("transform", translate(22, rectSize / 2))
  .attr("class", "legendText")
  .style("font-size", "10px")
  .attr("dy", ".35em")
  .text(function(d) { return d.label; });


///////////////////////////////////
// Legend mouseover and mouseout //
///////////////////////////////////

function selectLegend(opacity) {

  return function(d) {
    var chosen = d.range;

    wrapper.selectAll(".point")
      .filter(function(d) { return d.range != chosen; })
      .transition()
      .style("opacity", opacity);
  };

}

function selectShapeLegend(opacity) {

  return function(d) {
    var chosen = d.type;

    wrapper.selectAll(".point")
      .filter(function(d) { return d.type != chosen; })
      .transition()
      .style("opacity", opacity);
  };
}


//////////////////
// Legend click //
//////////////////

var chosenRange = new Set();
var showImputed = true,
    showGenotyped = true;

function clickShapeLegend(d) {

  d3.event.stopPropagation();

  var chosen = d.type;

  d3.selectAll(".legendShape")
    .style("opacity", function(d) {
      if (d.type !== chosen) return 0.3;
      else return 1;
    })
    .on("mouseover", null)
    .on("mouseout", null);

  if (chosen === "imputed") {
    showImputed = true;
    if (showGenotyped) showGenotyped = false;
    else {
      showGenotyped = true;
      d3.selectAll(".legendShape")
        .style("opacity", 1)
        .on("mouseover", selectShapeLegend(lowerOpacity))
        .on("mouseout", selectShapeLegend(1))
    }
  } else {
    showGenotyped = true;
    if (showImputed) showImputed = false;
    else {
      showImputed = true;
      d3.selectAll(".legendShape")
        .style("opacity", 1)
        .on("mouseover", selectShapeLegend(lowerOpacity))
        .on("mouseout", selectShapeLegend(1))
    }
  }

  // The imputed points
  wrapper.selectAll(".point.imputed")
    .style("opacity", 1)
    .style("visibility", function(d) {
      if ((showImputed) && ((chosenRange.size === 0) || (chosenRange.has(d.range)))) {
        return "visible";
      } else return "hidden";
    })
    .on("mouseover", function(d) {
      if ((showImputed) && ((chosenRange.size === 0) || (chosenRange.has(d.range)))) {
        return showTooltip.call(this, d);
      } else return null;
    })
    .on("mouseout", function(d) {
      if ((showImputed) && ((chosenRange.size === 0) || (chosenRange.has(d.range)))) {
        return removeTooltip.call(this, d);
      } else return null;
    })

  // The genotyped points
  wrapper.selectAll(".point.genotyped")
    .style("opacity", 1)
    .style("visibility", function(d) {
      if ((showGenotyped) && ((chosenRange.size === 0) || (chosenRange.has(d.range)))) {
        return "visible";
      } else return "hidden";
    })
    .on("mouseover", function(d) {
      if ((showGenotyped) && ((chosenRange.size === 0) || (chosenRange.has(d.range)))) {
        return showTooltip.call(this, d);
      } else return null;
    })
    .on("mouseout", function(d) {
      if ((showGenotyped) && ((chosenRange.size === 0) || (chosenRange.has(d.range)))) {
        return removeTooltip.call(this, d);
      } else return null;
    })

}

function clickLegend(d) {

  d3.event.stopPropagation();

  var chosen = d.range;

  if (chosenRange.has(chosen)) chosenRange.delete(chosen);
  else chosenRange.add(chosen);

  d3.selectAll(".legendSquare")
    .style("opacity", function(d) {
      if (chosenRange.has(d.range)) return 1;
      else return 0.3;
    })
    .on("mouseover", null)
    .on("mouseout", null);

  // The imputed points
  wrapper.selectAll(".point.imputed")
    .style("opacity", 1)
    .style("visibility", function(d) {
      if (!showImputed || (!chosenRange.has(d.range))) return "hidden";
      else return "visible";
    })
    .on("mouseover", function(d) {
      if (!showImputed || (!chosenRange.has(d.range))) return null;
      else return showTooltip.call(this, d);
    })
    .on("mouseout", function(d) {
      if (!showImputed || (!chosenRange.has(d.range))) return null;
      else return removeTooltip.call(this, d);
    });

  // The genotyped points
  wrapper.selectAll(".point.genotyped")
    .style("opacity", 1)
    .style("visibility", function(d) {
      if (!showGenotyped || (!chosenRange.has(d.range))) return "hidden";
      else return "visible";
    })
    .on("mouseover", function(d) {
      if (!showGenotyped || (!chosenRange.has(d.range))) return null;
      else return showTooltip.call(this, d);
    })
    .on("mouseout", function(d) {
      if (!showGenotyped || (!chosenRange.has(d.range))) return null;
      else return removeTooltip.call(this, d);
    });
}


/////////////////////
// The reset click //
/////////////////////

function resetClick() {

  chosenRange = new Set();
  showImputed = true;
  showGenotyped = true;

  d3.selectAll(".legendSquare")
    .style("opacity", 1)
    .on("mouseover", selectLegend(lowerOpacity))
    .on("mouseout", selectLegend(1));

  d3.selectAll(".legendShape")
    .style("opacity", 1)
    .on("mouseover", selectShapeLegend(lowerOpacity))
    .on("mouseout", selectShapeLegend(1));

  wrapper.selectAll(".point")
    .style("opacity", 1)
    .style("visibility", "visible");

  wrapper.selectAll(".point.imputed")
    .on("mouseover", showTooltip)
    .on("mouseout", removeTooltip);

  wrapper.selectAll(".point.genotyped")
    .on("mouseover", showTooltip)
    .on("mouseout", removeTooltip);

}

d3.select("body").on("click", resetClick);


