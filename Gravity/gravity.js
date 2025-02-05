﻿function zoom(direction) {
    if (direction > 0) {
        SCALE *= 1.1;
    }
    if (direction < 0) {
        SCALE /= 1.1;
    }
}

function panX(direction) {
    PAN_X += direction;
}

function panY(direction) {
    PAN_Y += direction;
}

function stopBall(op) {
    if (op == 1) {
        for (var i = 0; i < body.length; i++) {
            body[i].vx = 0;
            body[i].vy = 0;      
        }
    }
    if (op == 0) {
        for (var i = 0; i < body.length; i++) {
            if (body[i].still == true || body[i].grid == true) continue;
            body[i].x = Math.random()*GRID_WIDTH-GRID_WIDTH/2;
            body[i].y = Math.random()*GRID_HEIGHT-GRID_HEIGHT/2;
            body[i].vx = 0;
            body[i].vy = 0;      
        }
    }
    if (op == -1) {
        for (var i = 0; i < body.length; i++) {
            if (body[i].still == true) {
                body[i].x = Math.random()*GRID_WIDTH-GRID_WIDTH/2;
                body[i].y = Math.random()*GRID_HEIGHT-GRID_HEIGHT/2;
            }
        }
    }
}


function resetCanvas() {
    canvas.width = canvas.width;
    ctx.transform(SCALE, 0, 0, -SCALE, canvas.width / 2 + PAN_X, canvas.height / 2 + PAN_Y);
    ctx.lineWidth = 1 / SCALE;
}
resetCanvas();

//BODY PROPERTIES
var body = [];
for (var x = 0; x < N_G; x++) {
    for (var y = 0; y < N_G; y++) {
        body.push({
            color: "black",
            radius: 0.01,
            mass: 0.000001,
            grid: true,
            x: ((x + 0.5) / N_G - 0.5) * GRID_WIDTH,
            y: ((y + 0.5) / N_G - 0.5) * GRID_HEIGHT,
            vx: 0,
            vy: 0,
            ax: 0,
            ay: 0,
            fx: 0,
            fy: 0
        });
    }
}
body.push({
    color: "red",
    radius: 3,
    mass: 1,
    x: Math.random()*GRID_WIDTH-GRID_WIDTH/2,
    y: Math.random()*GRID_HEIGHT-GRID_HEIGHT/2,
    vx: 0,
    vy: 0,
    ax: 0,
    ay: 0,
    fx: 0,
    fy: 0      
})
function removeBody() {
    for (var i = 0; i < body.length; i++) {
        if(body[i].still == true && i == body.length - 1) {
            body.pop()
        }
    }
}
function placeBodies() {
    body.push({
        color: "yellow",
        mass: Math.random()*10*SUN_MASS,
        radius: 5,
        still: true,
        x: Math.random()*GRID_WIDTH-GRID_WIDTH/2,
        y: Math.random()*GRID_HEIGHT-GRID_HEIGHT/2,
        vx: 0,
        vy: 0,
        ax: 0,
        ay: 0,
        fx: 0,
        fy: 0
    });
}

//The gravitational force between exactly two objects.
function force(receiver, actor) {
    var dx = receiver.x - actor.x;
    var dy = receiver.y - actor.y;
    var rSquared = Math.pow(dx, 2) + Math.pow(dy, 2);
    var F = -(receiver.mass * actor.mass * G) / rSquared;
    var fx = F * dx / Math.sqrt(rSquared);
    var fy = F * dy / Math.sqrt(rSquared);
    return {
        x: fx,
        y: fy
    };
}

//The total force on every object.
function physics() {
    for (var i = 0; i < body.length; i++) {
        body[i].fx = 0;
        body[i].fy = 0;
        for (var j = 0; j < body.length; j++) {
            if (i == j) continue;
            var F = force(body[i], body[j]);
            body[i].fx += F.x;
            body[i].fy += F.y;
        }
    }
}


function moveObjects() {
    for (var i = 0; i < body.length; i++) {
        if (body[i].grid == true || body[i].still == true) continue;
        body[i].ax = body[i].fx / body[i].mass;
        body[i].ay = body[i].fy / body[i].mass;
        body[i].vx = body[i].vx + body[i].ax * dt;
        body[i].vy = body[i].vy + body[i].ay * dt;
        body[i].x = body[i].x + body[i].vx * dt;
        body[i].y = body[i].y + body[i].vy * dt;
    }
}

function drawField(gridPoint) {
    var F = Math.sqrt(Math.pow(gridPoint.fx, 2) + Math.pow(gridPoint.fy, 2));
    var lineLength = Math.sqrt((Math.log(F)/Math.log(10)+9)*F/1e-9);
    if (lineLength/SCALE > 1e11) {
        lineLength = 1e11*SCALE;
    }
    if (lineLength/SCALE < 0.5e10) {
        lineLength = 0.5e10*SCALE;
    }
    var dx = gridPoint.fx / F * lineLength/SCALE;
    var dy = gridPoint.fy / F * lineLength/SCALE;
    ctx.beginPath();
    ctx.moveTo(gridPoint.x - dx/2, gridPoint.y - dy/2);
    ctx.lineTo(gridPoint.x + dx/2, gridPoint.y + dy/2);
    ctx.strokeStyle = "black";
    ctx.stroke();
    ctx.beginPath();
    ctx.arc(gridPoint.x + dx/2, gridPoint.y + dy/2, 2/SCALE, 0, 2 * PI);
    ctx.stroke();
    //console.log(F, lineLength/SCALE)
  }

function drawObjects() {
    for (var i = 0; i < body.length; i++) {
        if (body[i].grid == true) {
            drawField(body[i]);
        } 
        if (body[i].still == true) {
            body[i].radius = body[i].mass/SUN_MASS;
        }
        ctx.beginPath();
        ctx.arc(body[i].x, body[i].y, body[i].radius / SCALE, 0, 2 * PI);
        ctx.strokeStyle = body[i].color;
        ctx.stroke();
        ctx.fillStyle = body[i].color;
        ctx.fill();
    }
}

function allOfExistence() {
    physics();
    moveObjects();
    resetCanvas();
    drawObjects();
}
setInterval(allOfExistence, SCALE_TIME * 1000 * dt);