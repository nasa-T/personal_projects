
function factorial(n) {
    var factorial = 1;
    for (var i = 1; i <= n; i++) {
        factorial *= i;
    }
    return factorial;
}

function derivative(f) {
    const dx = 0.1;
    return x => (f(x + dx) - f(x - dx)) / (2 * dx);
}

function f(x) {
    return x * Math.cos(x);
}

function nDerive(f, n) {
    var fPrime = f;
    for (var i = 0; i < n; i++) {
        fPrime = derivative(fPrime);
    }
    return fPrime;
}

function tSeries(x, n) {
    var sum = 0;
    for (var i = 0; i < n + 1; i++) {
        sum += nDerive(f, i)(SHIFT) / factorial(i) * (x - SHIFT) ** i;
    }
    if (sum > 1000) {
        sum = 1000;
    }
    return sum;
}

function resetCanvas() {
    canvas.width = canvas.width;
    ctx.transform(SCALE, 0, 0, -SCALE, canvas.width / 2, canvas.width / 2)
    ctx.lineWidth = 1 / SCALE;
}
resetCanvas();

function resetCanvas2() {
    canvas2.width = canvas2.width;
    ctx2.transform(SCALE, 0, 0, -SCALE, canvas2.width / 2, canvas2.width / 2)
    ctx2.lineWidth = 1 / SCALE;
    document.getElementById("degree").innerHTML = "degree: " + DEGREE;
    ctx2.beginPath();
    ctx2.moveTo(x, tSeries(x, DEGREE));
    ctx2.strokeStyle = "red";
}
resetCanvas2();

ctx.beginPath();
ctx.moveTo(-canvas.width / 2, 0);
ctx.lineTo(canvas.width / 2, 0);
ctx.stroke();
ctx.beginPath();
ctx.moveTo(0, -canvas.width / 2);
ctx.lineTo(0, canvas.width / 2);
ctx.stroke();

ctx.strokeStyle = "black";
ctx.fillStyle = "black";

function drawPlot1() {
    ctx.beginPath();
    ctx.arc(x, f(x), RADIUS / SCALE, 0, 2 * PI);
    ctx.stroke();
    //ctx.fill();
}

ctx2.beginPath();
ctx2.moveTo(x, tSeries(x, DEGREE));
var r = 255 * Math.random();
var g = 255 * Math.random();
var b = 255 * Math.random();
ctx2.strokeStyle = "rgb(" + r + ", " + g + ", " + b + ")";
ctx2.strokeStyle = "red";

function drawPlot2() {
    ctx2.lineTo(x, tSeries(x, DEGREE));
    ctx2.stroke();
}

for (x = DOMAIN_MIN; x < DOMAIN_MAX;) {
    drawPlot1();
    x += DOMAIN_LENGTH / N_PLOTS;
}

function allOfExistence() {
    x += DOMAIN_LENGTH / N_PLOTS;
    if (x <= DOMAIN_MAX) {
        tSeries();
        drawPlot2();
    } else if (DEGREE < MAX_DEGREE) {
        DEGREE += DEGREE_CHANGE;
        x = DOMAIN_MIN;
        resetCanvas2();
    }
}

setInterval(allOfExistence, 1000 * dt);
