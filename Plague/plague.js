var body = [];
var body2 = [];
for (var i = 0; i < N_B; i++) {
    var vx = Math.random() * 2 * SPEED - SPEED;
    var vy = (Math.round(Math.random()) * 2 - 1) * Math.sqrt(Math.pow(SPEED, 2) - Math.pow(vx, 2));
    body[i] = {
        color: "black",
        radius: 0.5,
        x: Math.random() * LENGTH - LENGTH / 2,
        y: Math.random() * HEIGHT - HEIGHT / 2,
        vx: vx,
        vy: vy,
        infected: false,
        dead: false,
        recovered: false,
        timeInfected: 0
    };
}
for (var i = 0; i < N_I; i++) {
    var vx = Math.random() * 2 * SPEED - SPEED;
    var vy = (Math.round(Math.random()) * 2 - 1) * Math.sqrt(Math.pow(SPEED, 2) - Math.pow(vx, 2));
    body[i] = {
        color: "red",
        radius: 0.5,
        x: Math.random() * LENGTH - LENGTH / 2,
        y: Math.random() * HEIGHT - HEIGHT / 2,
        vx: vx,
        vy: vy,
        infected: true,
        dead: false,
        recovered: false,
        timeInfected: 0
    };
}

for (var i = 0; i < N_B; i++) {
    var vx = Math.random() * 2 * SPEED - SPEED;
    var vy = (Math.round(Math.random()) * 2 - 1) * Math.sqrt(Math.pow(SPEED, 2) - Math.pow(vx, 2));
    body2[i] = {
        color: "black",
        radius: 0.5,
        x: Math.random() * LENGTH - LENGTH / 2,
        y: Math.random() * HEIGHT - HEIGHT / 2,
        vx: vx,
        vy: vy,
        infected: false,
        dead: false,
        recovered: false,
        timeInfected: 0
    };
}
for (var i = 0; i < N_I; i++) {
    var vx = Math.random() * 2 * SPEED - SPEED;
    var vy = (Math.round(Math.random()) * 2 - 1) * Math.sqrt(Math.pow(SPEED, 2) - Math.pow(vx, 2));
    body2[i] = {
        color: "red",
        radius: 0.5,
        x: Math.random() * LENGTH - LENGTH / 2,
        y: Math.random() * HEIGHT - HEIGHT / 2,
        vx: vx,
        vy: vy,
        infected: true,
        dead: false,
        recovered: false,
        timeInfected: 0
    };
}
var SCALE = canvas.width / (LENGTH);


function resetCanvas() {
    canvas.width = canvas.width;
    ctx.transform(SCALE, 0, 0, -SCALE, canvas.width / 2, canvas.height / 2);
    ctx.lineWidth = 1 / SCALE;
}

function resetCanvas2() {
    canvas2.width = canvas2.width;
    ctx2.transform(SCALE, 0, 0, -SCALE, canvas2.width / 2, canvas2.height / 2);
    ctx2.lineWidth = 1 / SCALE;
}
resetCanvas();
resetCanvas2();

function moveObjects() {
    for (var i = 0; i < body.length; i++) {
        if (body[i].dead == true) continue;
        body[i].x = body[i].x + body[i].vx * dt;
        body[i].y = body[i].y + body[i].vy * dt;
        if (body[i].x > LENGTH / 2) {
            body[i].x = LENGTH / 2;
            body[i].vx *= -1;
        }
        if (body[i].x < -LENGTH / 2) {
            body[i].x = -LENGTH / 2;
            body[i].vx *= -1;
        }
        if (body[i].y > HEIGHT / 2) {
            body[i].y = HEIGHT / 2;
            body[i].vy *= -1;
        }
        if (body[i].y < -HEIGHT / 2) {
            body[i].y = -HEIGHT / 2;
            body[i].vy *= -1;
        }
    }
}

function moveObjects2() {
    for (var i = 0; i < body2.length; i++) {
        if (body2[i].dead == true) continue;
        body2[i].x = body2[i].x + body2[i].vx * dt;
        body2[i].y = body2[i].y + body2[i].vy * dt;
        if (body2[i].x > LENGTH / 2) {
            body2[i].x = LENGTH / 2;
            body2[i].vx *= -1;
        }
        if (body2[i].x < -LENGTH / 2) {
            body2[i].x = -LENGTH / 2;
            body2[i].vx *= -1;
        }
        if (body2[i].y > HEIGHT / 2) {
            body2[i].y = HEIGHT / 2;
            body2[i].vy *= -1;
        }
        if (body2[i].y < -HEIGHT / 2) {
            body2[i].y = -HEIGHT / 2;
            body2[i].vy *= -1;
        }
        if (i % 2 == 0) {
            body2[i].vx = 0;
            body2[i].vy = 0;
        }
    }
}

function counting() {
    var count = 0;
    var dead = 0;
    var recovered = 0;
    var never = 0;
    for (var i = 0; i < body.length; i++) {
        if (body[i].infected == true) {
            count++;
        }
        if (body[i].dead == true) {
            dead++;
        }
        if (body[i].recovered == true) {
            recovered++;
        }
        if (body[i].infected == false && body[i].recovered == false && body[i].dead == false) {
            never++;
        }
    }
    displayText();
    return {
        infected: count,
        dead: dead,
        recovered: recovered,
        never: never,
        total: dead + recovered + count
    }

    function displayText() {
        document.getElementById("infected").innerHTML = count / body.length * 100 + "%";
        document.getElementById("n_infected").innerHTML = count + "of" + N_B;
        document.getElementById("dead").innerHTML = dead / body.length * 100 + "%"
        document.getElementById("n_dead").innerHTML = dead + "of" + N_B;
        document.getElementById("recovered").innerHTML = recovered / (recovered + dead) * 100 + "%"
        document.getElementById("n_recovered").innerHTML = recovered + "of" + (recovered + dead);
        document.getElementById("neverInf").innerHTML = never / body.length * 100 + "%"
        document.getElementById("n_neverInf").innerHTML = never + "of" + N_B;
        document.getElementById("total").innerHTML = (dead + recovered + count) / body.length * 100 + "%";
        document.getElementById("n_total").innerHTML = (dead + recovered + count) + "of" + N_B;
    }
}

function counting2() {
    var count = 0;
    var dead = 0;
    var recovered = 0;
    var never = 0;
    for (var i = 0; i < body2.length; i++) {
        if (body2[i].infected == true) {
            count++;
        }
        if (body2[i].dead == true) {
            dead++;
        }
        if (body2[i].recovered == true) {
            recovered++;
        }
        if (body2[i].infected == false && body2[i].recovered == false && body2[i].dead == false) {
            never++;
        }
    }
    displayText2();
    return {
        infected: count,
        dead: dead,
        recovered: recovered,
        never: never,
        total: dead + recovered + count
    }

    function displayText2() {
        document.getElementById("infected2").innerHTML = count / body2.length * 100 + "%";
        document.getElementById("n_infected2").innerHTML = count + "of" + N_B;
        document.getElementById("dead2").innerHTML = dead / body2.length * 100 + "%"
        document.getElementById("n_dead2").innerHTML = dead + "of" + N_B;
        document.getElementById("recovered2").innerHTML = recovered / (recovered + dead) * 100 + "%"
        document.getElementById("n_recovered2").innerHTML = recovered + "of" + (recovered + dead);
        document.getElementById("neverInf2").innerHTML = never / body2.length * 100 + "%"
        document.getElementById("n_neverInf2").innerHTML = never + "of" + N_B;
        document.getElementById("total2").innerHTML = (dead + recovered + count) / body2.length * 100 + "%";
        document.getElementById("n_total2").innerHTML = (dead + recovered + count) + "of" + N_B;
    }
}

function infection() {
    for (var i = 0; i < body.length; i++) {
        if (body[i].dead == true) continue;
        for (var j = 0; j < body.length; j++) {
            if (i == j) continue;
            if (Math.abs(body[i].x - body[j].x) < (body[i].radius + body[j].radius) + ACHOO && Math.abs(body[i].y - body[j].y) < (body[i].radius + body[j].radius) + ACHOO && body[j].infected == true) {
                body[i].infected = true;
                body[i].color = "red";
            }
        }
        if (body[i].infected == true) {
            body[i].timeInfected += 1;
        }
    }
}

function infection2() {
    for (var i = 0; i < body2.length; i++) {
        if (body2[i].dead == true) continue;
        for (var j = 0; j < body2.length; j++) {
            if (i == j) continue;
            if (Math.abs(body2[i].x - body2[j].x) < (body2[i].radius + body2[j].radius) + ACHOO && Math.abs(body2[i].y - body2[j].y) < (body2[i].radius + body2[j].radius) + ACHOO && body2[j].infected == true) {
                body2[i].infected = true;
                body2[i].color = "red";
            }
        }
        if (body2[i].infected == true) {
            body2[i].timeInfected += 1;
        }
    }
}


function recovery() {
    for (var i = 0; i < body.length; i++) {
        if (body[i].dead == false && body[i].infected == true) {
            if (body[i].timeInfected > TIME_FATE) {
                var chance = Math.random() * 100;
                if (chance < FATE_DECIDER) {
                    body[i].recovered = true;
                    body[i].infected = false;
                }
            }
            if (body[i].recovered == true) {
                body[i].infected = false;
            } //reinfection is not allowed
            if (body[i].recovered == true) {
                body[i].color = "green";
            }
        }
    }
}

function recovery2() {
    for (var i = 0; i < body2.length; i++) {
        if (body2[i].dead == false && body2[i].infected == true) {
            if (body2[i].timeInfected > TIME_FATE) {
                var chance = Math.random() * 100;
                if (chance < FATE_DECIDER) {
                    body2[i].recovered = true;
                    body2[i].infected = false;
                }
            }
            if (body2[i].recovered == true) {
                body2[i].infected = false;
            } //reinfection is not allowed
            if (body2[i].recovered == true) {
                body2[i].color = "green";
            }
        }
    }
}


function death() {
    for (var i = 0; i < body.length; i++) {
        if (body[i].timeInfected > TIME_FATE && body[i].infected == true) {
            var chance = Math.random() * 100;
            if (chance > FATE_DECIDER) {
                body[i].dead = true;
            }
        }
        if (body[i].dead == true) {
            body[i].infected = false;
        }
    }
}

function death2() {
    for (var i = 0; i < body2.length; i++) {
        if (body2[i].timeInfected > TIME_FATE && body2[i].infected == true) {
            var chance = Math.random() * 100;
            if (chance > FATE_DECIDER) {
                body2[i].dead = true;
            }
        }
        if (body2[i].dead == true) {
            body2[i].infected = false;
        }
    }
}

function drawObjects() {
    for (var i = 0; i < body.length; i++) {
        if (body[i].dead == false) {
            ctx.beginPath();
            ctx.arc(body[i].x, body[i].y, body[i].radius, 0, 2 * PI);
            ctx.strokeStyle = body[i].color;
            ctx.stroke();
            ctx.fillStyle = body[i].color;
            ctx.fill();
        }
    }
}
function drawObjects2() {
    for (var i = 0; i < body2.length; i++) {
        if (body2[i].dead == false) {
            ctx2.beginPath();
            ctx2.arc(body2[i].x, body2[i].y, body2[i].radius, 0, 2 * PI);
            ctx2.strokeStyle = body2[i].color;
            ctx2.stroke();
            ctx2.fillStyle = body2[i].color;
            ctx2.fill();
        }
    }
}

var GRAPH_SCALE = graph.height / N_B;

function resetGraph() {
    graph.width = graph.width;
    con.transform(GRAPH_SCALE / 6, 0, 0, -GRAPH_SCALE, 0, graph.height);
    con.lineWidth = 1 / GRAPH_SCALE;
}
function resetGraph2() {
    graph2.width = graph2.width;
    con2.transform(GRAPH_SCALE / 6, 0, 0, -GRAPH_SCALE, 0, graph2.height);
    con2.lineWidth = 1 / GRAPH_SCALE;
}

var total = [];
var infected = [];
var mortality = [];
var recover = [];
var untouched = [];
function drawGraph() {
    var C = counting();
    total.unshift({
        x: t,
        y: C.total,
        color: "purple"
    });
    infected.unshift({
        x: t,
        y: C.infected,
        color: "red"
    });
    mortality.unshift({
        x: t,
        y: C.dead,
        color: "white"
    });
    recover.unshift({
        x: t,
        y: C.recovered,
        color: "green"
    });
    untouched.unshift({
        x: t,
        y: C.never,
        color: "black"
    });
    var i = 0
    con.beginPath();
    con.arc(total[i].x, total[i].y, PLOT_SIZE / GRAPH_SCALE, 0, 2 * PI);
    con.strokeStyle = total[i].color;
    con.stroke();
    con.fillStyle = total[i].color;
    con.fill();
    con.beginPath();
    con.arc(infected[i].x, infected[i].y, PLOT_SIZE / GRAPH_SCALE, 0, 2 * PI);
    con.strokeStyle = infected[i].color;
    con.stroke();
    con.fillStyle = infected[i].color;
    con.fill();
    con.beginPath();
    con.arc(mortality[i].x, mortality[i].y, PLOT_SIZE / GRAPH_SCALE, 0, 2 * PI);
    con.strokeStyle = mortality[i].color;
    con.stroke();
    con.fillStyle = mortality[i].color;
    con.fill();
    con.beginPath();
    con.arc(recover[i].x, recover[i].y, PLOT_SIZE / GRAPH_SCALE, 0, 2 * PI);
    con.strokeStyle = recover[i].color;
    con.stroke();
    con.fillStyle = recover[i].color;
    con.fill();
    con.beginPath();
    con.arc(untouched[i].x, untouched[i].y, PLOT_SIZE / GRAPH_SCALE, 0, 2 * PI);
    con.strokeStyle = untouched[i].color;
    con.stroke();
    con.fillStyle = untouched[i].color;
    con.fill();
}
resetGraph();

var total2 = [];
var infected2 = [];
var mortality2 = [];
var recover2 = [];
var untouched2 = [];
function drawGraph2() {
    var C2 = counting2();
    total2.unshift({
        x: t,
        y: C2.total,
        color: "purple"
    });
    infected2.unshift({
        x: t,
        y: C2.infected,
        color: "red"
    });
    mortality2.unshift({
        x: t,
        y: C2.dead,
        color: "white"
    });
    recover2.unshift({
        x: t,
        y: C2.recovered,
        color: "green"
    });
    untouched2.unshift({
        x: t,
        y: C2.never,
        color: "black"
    });
    var i = 0
    con2.beginPath();
    con2.arc(total2[i].x, total2[i].y, PLOT_SIZE / GRAPH_SCALE, 0, 2 * PI);
    con2.strokeStyle = total2[i].color;
    con2.stroke();
    con2.fillStyle = total2[i].color;
    con2.fill();
    con2.beginPath();
    con2.arc(infected2[i].x, infected2[i].y, PLOT_SIZE / GRAPH_SCALE, 0, 2 * PI);
    con2.strokeStyle = infected2[i].color;
    con2.stroke();
    con2.fillStyle = infected2[i].color;
    con2.fill();
    con2.beginPath();
    con2.arc(mortality2[i].x, mortality2[i].y, PLOT_SIZE / GRAPH_SCALE, 0, 2 * PI);
    con2.strokeStyle = mortality2[i].color;
    con2.stroke();
    con2.fillStyle = mortality2[i].color;
    con2.fill();
    con2.beginPath();
    con2.arc(recover2[i].x, recover2[i].y, PLOT_SIZE / GRAPH_SCALE, 0, 2 * PI);
    con2.strokeStyle = recover2[i].color;
    con2.stroke();
    con2.fillStyle = recover2[i].color;
    con2.fill();
    con2.beginPath();
    con2.arc(untouched2[i].x, untouched2[i].y, PLOT_SIZE / GRAPH_SCALE, 0, 2 * PI);
    con2.strokeStyle = untouched2[i].color;
    con2.stroke();
    con2.fillStyle = untouched2[i].color;
    con2.fill();
}
resetGraph2();

function allOfExistence() {
    t += 1;
    moveObjects();
    moveObjects2();
    infection();
    infection2();
    recovery();
    recovery2();
    death();
    death2();
    if (t % 5 == 0) {
        resetCanvas();
        resetCanvas2();
        drawObjects();
        drawObjects2();
        drawGraph();
        drawGraph2();
    }
}
setInterval(allOfExistence, SCALE_TIME * 1000 * dt);
