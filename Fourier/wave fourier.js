var wave = document.getElementById("canvas");
var ctx = wave.getContext("2d");
var bar = document.getElementById("graph");
var con = bar.getContext("2d");

console.log(WAVE.length)

const PI = Math.PI;
const dt = 0.01;
const DOMAIN_MIN = 0;
const DOMAIN_MAX = WAVE.length;
const WAVE_MAX = 100; // max range of wave
const SCALE = wave.height / WAVE_MAX;
const N_MAX = 100;
const G_SCALE = bar.width / (N_MAX);
const N_INTERVALS = WAVE.length - 1;
var xMin = 0;
var xMax = wave.width;
var T = (xMax - xMin) / 2;


function resetCanvas() {
	wave.width = wave.width;
	ctx.transform(1, 0, 0, -SCALE, -xMin, 0);
	ctx.lineWidth = 1 / (SCALE);
	bar.width = bar.width;
	con.transform(G_SCALE, 0, 0, -bar.height / 10, 0, bar.height)
	con.lineWidth = 1 / (bar.height / 1);
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
	for (var x = xMin; x < xMax; x++) {
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

let a_n = [];

let b_n = [];

function coefficiate(f) {
	var total = 0;
	for (var i = 0; i < N_MAX; i++) {
		a_n[i] = a(f, i);
		b_n[i] = b(f, i);
		total += Math.abs(a_n[i] + b_n[i]);
		//console.log((Math.sqrt(a_n[i] ** 2 + b_n[i] ** 2)).toFixed(3));
	}

	//console.log(total);
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
		var total = Math.abs(a_n[i] + b_n[i]);
		con.beginPath();
		con.fillRect(i, 0, 0.7, total);
		con.stroke();
	}
}

function drawAxes() {
	con.beginPath();
	con.moveTo(0, 0);
	con.lineTo(bar.width, 0);
	con.stroke();
}
drawAxes();

ctx.beginPath();
ctx.moveTo(DOMAIN_MIN, f(DOMAIN_MIN));
function drawGraph(x) {
	ctx.strokeStyle = "green";
	ctx.lineTo(x, f(x));
	ctx.stroke();
}

// build a wave
// var wave = [];
// for (var i = DOMAIN_MIN; i <= DOMAIN_MAX; i += (X_MAX - X_MIN) / N_INTERVALS) {
// 	wave.push(f(i));
// 	console.log("wave:", wave);
// }

// draw bar graph
coefficiate(f);
drawBars();

// draw function graphs


setInterval(
	() => {
		const SPEED = 1000;
		if (xMax < WAVE.length) {
			resetCanvas();
			for (var i = xMin; i < xMax; i++) {
				drawGraph(i);
			}
			drawAxes();
			coefficiate(f);
			drawBars();
			if (xMax + SPEED < WAVE.length) {
				xMin += SPEED;
				xMax += SPEED;
			} else {
				xMin += WAVE.length - xMax;
				xMax += WAVE.length - xMax;
			}
		}
	},
	1000 * dt
)
