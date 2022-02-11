//BODY PROPERTIES

var body = [];

for (var i = 0; i < N_E; i++) {
    body[i] = {
        color: "blue",
        mass: MASS_P,
        charge: e,
        radius: 1,
        wire: true,
        x: LENGTH * Math.random() - (LENGTH * 0.5),
        y: HEIGHT * Math.random() - (HEIGHT * 0.5),
        z: 0,
        vx: 0,
        vy: 0,
        vz: 0,
        ax: 0,
        ay: 0,
        az: 0,
        fx: 0,
        fy: 0,
        fz: 0
    }
}
for (var x = 0; x < N_G; x++) {
    for (var y = 0; y < N_G; y++) {
        body.push({
            color: "black",
            radius: 0.1,
            mass: Infinity,
            charge: e,
            grid: true,
            x: ((x + 0.5) / N_G - 0.5) * GRID_WIDTH,
            y: ((y + 0.5) / N_G - 0.5) * GRID_HEIGHT,
            z: 0,
            vx: 0,
            vy: 0,
            ax: 0,
            ay: 0,
            fx: 0,
            fy: 0
       });
   }
}

// for (var i = 0; i < N_P; i++) {
//     body.push({
//         color: "blue",
//         mass: MASS_P,
//         charge: e,
//         radius: 1,
//         still: false,
//         circle: true,
//         x: LENGTH * Math.random() - (LENGTH * 0.5),
//         y: HEIGHT * Math.random() - (HEIGHT * 0.5),
//         z: 0,
//         vx: 0,
//         vy: 0,
//         vz: 0,
//         ax: 0,
//         ay: 0,
//         az: 0,
//         fx: 0,
//         fy: 0,
//         fz: 0
//     })
// }

// body.push({
//     color: "black",
//     mass: MASS_E,
//     charge: e,
//     radius: 5,
//     still: false,
//     x: GRID_WIDTH * Math.random() - GRID_WIDTH/2,
//     y: GRID_HEIGHT * Math.random() - GRID_WIDTH/2,
//     vx: 0,
//     vy: 0,
//     ax: 0,
//     ay: 0,
//     fx: 0,
//     fy: 0
// })



var SCALE = canvas.width / (GRID_WIDTH);

function resetCanvas() {
    canvas.width = canvas.width;
    ctx.transform(SCALE, 0, 0, -SCALE, canvas.width / 2, canvas.height / 2);
    ctx.lineWidth = 1 / SCALE;
}
resetCanvas();

//The gravitational force between exactly two objects.
function force(receiver, actor) {
    if (actor.grid == true) {
        return {
            x: 0,
            y: 0,
            z: 0
        }
    }
    var dx = receiver.x - actor.x;
    var dy = receiver.y - actor.y;
    var dz = receiver.z - actor.z;
    var rSquared = Math.pow(dx, 2) + Math.pow(dy, 2) + Math.pow(dz, 2);
    if (rSquared < 0.000001) rSquared = 0.000001;
    var F = (K * receiver.charge * actor.charge) / rSquared;
    var fx = F * dx / Math.sqrt(rSquared);
    var fy = F * dy / Math.sqrt(rSquared);
    var fz = F * dz / Math.sqrt(rSquared);
    return {
        x: fx,
        y: fy,
        z: fz
    };
}

//The total force on every object.
function physics() {
    for (var i = 0; i < body.length; i++) {
        body[i].fx = 0;
        body[i].fy = 0;
        body[i].fz = 0;
        for (var j = 0; j < body.length; j++) {
            if (i == j) continue;
            var F = force(body[i], body[j]);
            body[i].fx += F.x;
            body[i].fy += F.y;
            body[i].fz += F.z;
        }
        
        if (body[i].grid == false) {
            if (body[i].fx > TOO_MUCH) body[i].fx = TOO_MUCH;
            if (body[i].fx < -TOO_MUCH) body[i].fx = -TOO_MUCH;
        }
    }
}


function moveObjects() {
    for (var i = 0; i < body.length; i++) {
        if (body[i].grid == true) {
            continue;
        }
        body[i].fz = 0;
        body[i].vz = 0;
        body[i].z = 0;
        body[i].ax = body[i].fx / body[i].mass;
        body[i].ay = body[i].fy / body[i].mass;
        body[i].az = body[i].fz / body[i].mass;
        body[i].vx = body[i].vx + body[i].ax * dt;
        body[i].vy = body[i].vy + body[i].ay * dt;
        body[i].vz = body[i].vz + body[i].az * dt;
        body[i].x = body[i].x + body[i].vx * dt;
        body[i].y = body[i].y + body[i].vy * dt;
        //body[i].z = body[i].z + body[i].vz * dt;
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
        if (body[i].z > DEPTH / 2) {
            body[i].z = DEPTH / 2;
            body[i].vz *= -1;
        }
        if (body[i].z < -DEPTH / 2) {
            body[i].z = -DEPTH / 2;
            body[i].vz *= -1;
        }
        body[i].vx *= FRICTION;
        body[i].vy *= FRICTION;
        body[i].vz *= FRICTION;
        if (body[i].wire == true) {
            body[i].y = 0;
            body[i].z = 0;
        }
        if (body[i].circle == true) {
            
        }
        if (body[i].still == true) {
            body[i].vy = 0;
            body[i].vx = 0;
            body[i].vz = 0;
        }
    }
}

function drawField(gridPoint) {
    var F = Math.sqrt(Math.pow(gridPoint.fx, 2) + Math.pow(gridPoint.fy, 2));
    // console.log(F)
    var lineLength = (Math.log(F) / Math.log(10) + 31) * 4;
    var dx = gridPoint.fx / F * lineLength / SCALE;
    var dy = gridPoint.fy / F * lineLength / SCALE;
    ctx.beginPath();
    ctx.moveTo(gridPoint.x - dx, gridPoint.y - dy);
    ctx.lineTo(gridPoint.x + dx, gridPoint.y + dy);
    ctx.strokeStyle = "black";
    ctx.stroke();
    ctx.beginPath();
    ctx.arc(gridPoint.x + dx, gridPoint.y + dy, 2 / SCALE, 0, 2 * PI);
    ctx.stroke();
}

function drawObjects() {
    for (var i = 0; i < body.length; i++) {
        if (body[i].grid == true) {
            drawField(body[i]);
        }
        ctx.beginPath();
        ctx.arc(body[i].x, body[i].y, (body[i].radius / SCALE), 0, 2 * PI);
        ctx.strokeStyle = body[i].color;
        ctx.stroke();
        ctx.fillStyle = body[i].color;
        ctx.fill();
    }
    ctx.beginPath();
    ctx.moveTo(-LENGTH / 2, 0);
    ctx.lineTo(LENGTH / 2, 0);
    ctx.strokeStyle = "green";
    ctx.stroke();
}

function allOfExistence() {
    physics();
    moveObjects();
    resetCanvas();
    drawObjects();
}
// allOfExistence()

setInterval(allOfExistence, SCALE_TIME * 1000 * dt);
