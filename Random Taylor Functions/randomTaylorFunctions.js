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
    dt = 0.1 / DEGREE;
    if (dt < 0.01) {
        dt = 0.01;
    }
}
resetCanvas2();

function drawAxes() {
    ctx.beginPath();
    ctx.moveTo(-canvas.width / 2, 0);
    ctx.lineTo(canvas.width / 2, 0);
    ctx.stroke();
    ctx.beginPath();
    ctx.moveTo(0, -canvas.width / 2);
    ctx.lineTo(0, canvas.width / 2);
    ctx.stroke();
}
drawAxes();
ctx.strokeStyle = "black";
ctx.fillStyle = "black";

function drawPlot1(x) {
    ctx.beginPath();
    ctx.arc(x, f(x), RADIUS / SCALE, 0, 2 * PI);
    ctx.stroke();
    //ctx.fill();
}

ctx2.beginPath();
ctx2.moveTo(DOMAIN_MIN, tSeries(DOMAIN_MIN, DEGREE));
var r = 255 * Math.random();
var g = 255 * Math.random();
var b = 255 * Math.random();
ctx2.strokeStyle = "rgb(" + r + ", " + g + ", " + b + ")";

function drawPlot2(x) {
    ctx2.lineTo(x, tSeries(x, DEGREE));
    ctx2.strokeStyle = "red";
    ctx2.stroke();
}

function toggleRun(button) {
    var toggleRun = document.getElementById("toggleRun");
    if (toggleRun.innerHTML == "Pause") {
        clearInterval(intervalId);
        toggleRun.innerHTML = "Restart";
    } else {
        clearInterval(intervalId);
        intervalId = setInterval(allOfExistence, 1000 * dt);
        toggleRun.innerHTML = "Pause";
    }
}

for (var x = DOMAIN_MIN; x < DOMAIN_MAX;) {
    drawPlot1(x);
    x += DOMAIN_LENGTH / N_PLOTS;
}

function allOfExistence() {
    SHIFT += DOMAIN_LENGTH / N_PLOTS;
    var fShift = f(SHIFT);
    while (
        SHIFT < DOMAIN_MAX &&
        (
            isNaN(fShift) ||
            fShift > DOMAIN_MAX * 100 ||
            fShift < DOMAIN_MIN * 100
        )
    ) {
        SHIFT += DOMAIN_LENGTH / N_PLOTS;
        fShift = f(SHIFT);
    }
    if (SHIFT <= DOMAIN_MAX) {
        resetCanvas2();
        var x = DOMAIN_MIN;
        while (x < DOMAIN_MAX) {
            drawPlot2(x);
            x += DOMAIN_LENGTH / N_PLOTS;
        }
        ctx2.beginPath();
        ctx2.arc(SHIFT, fShift, 2 / SCALE, 0, 2 * PI);
        ctx2.stroke();
        ctx2.fillStyle = "red";
        ctx2.fill();
    } else if (DEGREE < MAX_DEGREE) {
        SHIFT = DOMAIN_MIN;
        DEGREE += DEGREE_CHANGE;
    }
}

var intervalId = setInterval(allOfExistence, 400 * dt);
document.getElementById("toggleRun").innerHTML = "Pause";
