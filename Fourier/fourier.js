var canvas = document.getElementById("canvas");
var ctx = canvas.getContext("2d");
var canvas2 = document.getElementById("canvas2");
var ctx2 = canvas2.getContext("2d");
var graph = document.getElementById("graph");
var con = graph.getContext("2d");

console.log(WAVE.length)

const PI = Math.PI;
const dt = 0.01;
const DOMAIN_MIN = 0;
const DOMAIN_MAX = WAVE.length;
const SCALE = canvas.width / (DOMAIN_MAX - DOMAIN_MIN);
const N_MAX = 100;
const G_SCALE = canvas.width / (N_MAX + 2);
const N_INTERVALS = WAVE.length - 1;
var xMin = 0;
var xMax = WAVE.length;
const T = (xMax - xMin) * 2 / 2;


function resetCanvas() {
	canvas.width = canvas.width;
	ctx.transform(SCALE, 0, 0, -canvas.height / 200, 0, canvas.height / 2);
	ctx.lineWidth = 1 / (canvas.height / 200);
	canvas2.width = canvas2.width;
	ctx2.transform(SCALE, 0, 0, -canvas2.height / 200, 0, canvas2.height / 2);
	ctx2.lineWidth = 1 / (canvas2.height / 200);
	graph.width = graph.width;
	con.transform(G_SCALE, 0, 0, -graph.height / 1, 0, graph.height)
	con.lineWidth = 1 / G_SCALE;
}
resetCanvas();

function f(x) {
	return WAVE[x];
}

function integral(f, xMin = -PI, xMax = PI) {
	const dx = (xMax - xMin) / N_INTERVALS;
	var sum = 0;
	// for (var x = xMin + dx / 2; x < xMax; x += dx) {
	// 	sum += f(x) * dx;
	// }
	for (var x = xMin; x < xMax; x++){
		sum += f(x);
	}
	return sum;
}

function a(f, n) {
	if (n == 0) {
		return 1 / (2 * T) * integral(f, xMin, xMax);
	} else {
		const g = x => f(x) * Math.cos(n * PI / T * x);
		return 1 / T * integral(g, xMin, xMax);
	}
}

function b(f, n) {
	const g = x => f(x) * Math.sin(n * PI / T * x);
	return 1 / T * integral(g, xMin, xMax);
}
console.log(integral(f, xMin, xMax), f(0))

let a_n = [];

let b_n = [];

function coefficiate(f) {
	var total = 0;
	for (var i = 0; i <= N_MAX; i++) {
		a_n[i] = a(f, i);
		b_n[i] = b(f, i);
		total += a_n[i] ** 2 + b_n[i] ** 2;
		console.log(a_n[i].toFixed(3), b_n[i].toFixed(3));
	}
	
	console.log(total);
}

function fSeries(f, x) {
	var series = 0;
	for (var i = 0; i <= N_MAX; i++) {
		series += a_n[i] * Math.cos(i * PI / T * x) + b_n[i] * Math.sin(i * PI / T * x);
	}
	return series;
}

function drawBars() {
	for (var i = 0; i <= N_MAX; i++) {
		var total = Math.sqrt(a_n[i] ** 2 + b_n[i] ** 2);
		con.beginPath();
		con.fillRect(i, 0, 1, total);
		con.stroke();
	}
}

function drawAxes() {
	ctx.beginPath();
	ctx.moveTo(0, 0);
	ctx.lineTo(canvas.width / SCALE, 0);
	ctx.stroke();
	ctx.beginPath();
	ctx.moveTo(0, -canvas.height / SCALE);
	ctx.lineTo(0, canvas.height / SCALE);
	ctx.stroke();
}
drawAxes();

ctx.beginPath();
ctx2.beginPath();
ctx.moveTo(DOMAIN_MIN, f(DOMAIN_MIN));
//ctx2.moveTo(DOMAIN_MIN, fSeries(f, DOMAIN_MIN));
ctx.strokeStyle = "green";
ctx2.strokeStyle = "red";
function drawGraph(x) {
	ctx.lineTo(x, f(x));
	//ctx2.lineTo(x, fSeries(f, x));
	ctx.stroke();
	ctx2.stroke();
}

// var wave = [];


// for (var i = DOMAIN_MIN; i <= DOMAIN_MAX; i += (X_MAX - X_MIN) / N_INTERVALS) {
// 	wave.push(f(i));
// 	console.log("wave:", wave);
// }

// draw bar graph
coefficiate(f);
drawBars();

// draw function graphs
for (var i = DOMAIN_MIN; i < DOMAIN_MAX; i += 5) {
	drawGraph(i);
}

// setInterval(
// 	() => {
		// if (x <= DOMAIN_MAX) {
		// 	drawGraph(x);
		// 	x += (DOMAIN_MAX - DOMAIN_MIN) / N_INTERVALS;
		// }
// 	},
// 	1000 * dt
// )
