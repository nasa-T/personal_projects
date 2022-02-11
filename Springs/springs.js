var pNum = Math.floor(Math.random() * parameters.length);
//pNum = 2;
var p = parameters[pNum];
const SCALE = p.SCALE;

var point = [];

for (var i = 0; i < p.nJoints; i++) {
    point[i] = {
        x: 0,
        y: (p.start - i * p.springLength),
        vx: 0,
        vy: 0,
        mass: 1,
        k: p.k
    };

    if (p.sideways == true) {
        point[i] = {
            y: 0,
            x: -p.start + i * (p.start * 2) / (p.nJoints - 1),
            vx: 0,
            vy: 0,
            mass: 1,
            k: p.k
        };
    }

    if (p.attach == true) {
        point[i] = {
            x: 0,
            y: p.start - i * (p.start * 2) / (p.nJoints - 1),
            vx: 0,
            vy: 0,
            mass: 1,
            k: p.k
        };
    }

    if (i > 0) Object.assign(point[i], p.notTop);
    if (i == 1) Object.assign(point[i], p.second);
    if (i == 2) Object.assign(point[i], p.third);
    if (i == p.nJoints - 1) Object.assign(point[i], p.bottom);
    Object.assign(point[0], p.top);
}
console.log(pNum);

function resetCanvas() {
    canvas.width = canvas.width;
    ctx.transform(SCALE, 0, 0, -SCALE, canvas.width / 2, canvas.width / 2);
    ctx.lineWidth = 1 / SCALE;
}

(function resetCanvas2() {
    trail.width = trail.width;
    ctx2.transform(SCALE, 0, 0, -SCALE, trail.width / 2, trail.width / 2);
    ctx2.lineWidth = 1 / SCALE;
})();

function force(top, middle, bottom) {
    if (top == undefined) {
        var fGravity = -middle.mass * g;
        var dxb = bottom.x - middle.x;
        var dyb = bottom.y - middle.y;
        var dLb = Math.sqrt(dxb ** 2 + dyb ** 2);
        var stretchb = dLb - (bottom.mass * g + bottom.k * p.springLength) / bottom.k;
        var fSpringb = -bottom.k * stretchb;
        return {
            x: -fSpringb * (dxb / dLb),
            y: -fSpringb * (dyb / dLb) + fGravity
        }
    }
    if (bottom == undefined) {
        var fGravity = -middle.mass * g;
        var dx = middle.x - top.x;
        var dy = middle.y - top.y;
        var dL = Math.sqrt(dx ** 2 + dy ** 2);
        var stretch = dL - (middle.mass * g + middle.k * p.springLength) / middle.k;
        var fSpring = -middle.k * stretch;
        return {
            x: fSpring * (dx / dL),
            y: fSpring * (dy / dL) + fGravity
        }
    }
    var fGravity = -middle.mass * g;
    var dxm = middle.x - top.x;
    var dym = middle.y - top.y;
    var dLm = Math.sqrt(dxm ** 2 + dym ** 2);
    var stretchm = dLm - (middle.mass * g + middle.k * p.springLength) / middle.k;
    var fSpringm = -middle.k * stretchm;
    var dxb = bottom.x - middle.x;
    var dyb = bottom.y - middle.y;
    var dLb = Math.sqrt(dxb ** 2 + dyb ** 2);
    var stretchb = dLb - (bottom.mass * g + bottom.k * p.springLength) / middle.k;
    var fSpringb = -bottom.k * stretchb;
    return {
        x: fSpringm * (dxm / dLm) - fSpringb * (dxb / dLb),
        y: fSpringm * (dym / dLm) - fSpringb * (dyb / dLb) + fGravity
    }

}

function physics() {
    for (var i = 0; i < p.nJoints; i++) {
        var F = force(point[i - 1], point[i], point[i + 1]);
        point[i].fx = F.x;
        point[i].fy = F.y;
    }
}

var vx = p.shakeSpeed;

function movePoint() {
    for (var i = 0; i < p.nJoints; i++) {
        if (p.oneShake == false && t < p.timeToDrop) {
            point[i].vy *= p.friction;
            point[0].y = p.start;
            point[0].fx = 0;
            point[0].fy = 0;
            if (t > 500) {
                point[0].vx = vx;
                if (point[0].x > SHAKE_SIZE) {
                    point[0].x = SHAKE_SIZE;
                    vx = -p.shakeSpeed;
                }
                if (point[0].x < -SHAKE_SIZE) {
                    point[0].x = -SHAKE_SIZE;
                    vx = p.shakeSpeed;
                }
            }
        }
        if (p.oneShake == true && t < p.timeToDrop) {
            point[i].vy *= p.friction;
            point[0].y = p.start;
            point[0].fx = 0;
            point[0].fy = 0;
            if (t > 500) {
                point[0].vx = vx;
                if (point[0].x > SHAKE_SIZE) {
                    point[0].x = SHAKE_SIZE;
                    vx = -p.shakeSpeed;
                }
                if (point[0].x < -SHAKE_SIZE) {
                    point[0].x = -SHAKE_SIZE;
                    vx = p.shakeSpeed;
                }
                if (point[0].x < 0) {
                    point[0].vx = 0;
                }
            }
        }
        if (p.sideways == true) {
            point[0].x = -p.start;
            point[0].y = 0;
            point[p.nJoints - 1].x = p.start;
            point[p.nJoints - 1].vx = 0;
            point[i].y = 0;
            point[i].vy = 0;
        }
        point[i].ax = point[i].fx / point[i].mass;
        point[i].vx += point[i].ax * dt;
        point[i].x += point[i].vx * dt;
        point[i].ay = point[i].fy / point[i].mass;
        point[i].vy += point[i].ay * dt;
        point[i].y += point[i].vy * dt;
        if (i < p.nJoints - 1 && point[i].y < point[i + 1].y) {
            point[i].vy = point[i + 1].vy;
        }
        if (p.attach == true && i == p.nJoints - 1) {
            point[i].y = -p.start;
            point[i].x = 0;
        }
    }
}

function drawSprings() {
    for (var i = 0; i < p.nJoints; i++) {
        if (p.trace == true) {
            ctx2.beginPath();
            ctx2.arc(point[i].x, point[i].y, 1 / SCALE, 0, 2 * PI);
            ctx2.fill();
            ctx2.stroke();
        }
        if (i == p.nJoints - 1) continue;
        var dx = (point[i + 1].x - point[i].x) / p.nCoils;
        var dy = (point[i + 1].y - point[i].y) / p.nCoils;
        var width = p.width / Math.pow(Math.abs(point[0].y - point[p.nJoints - 1].y), 1 / 10);
        if (width > p.width) {
            width = p.width;
        }
        ctx.beginPath();
        ctx.moveTo(point[i].x, point[i].y);
        for (var j = 0; j < p.nCoils; j++) {
            if (p.sideways == false) {
                ctx.lineTo(point[i].x + width / 2 + dx * j, point[i].y + dy * j);
                if (j < p.nCoils - 1) {
                    ctx.lineTo(point[i].x - width / 2 + dx * j, point[i].y + dy * (j + 1));
                }
                if (j == p.nCoils - 1) {
                    ctx.lineTo(point[i].x - width / 2 + dx * j, point[i].y + dy * (j + 1));
                    ctx.lineTo(point[i + 1].x, point[i + 1].y);
                }
                ctx.stroke();
            }

            if (p.sideways == true) {
                ctx.lineTo(point[i].x + dx * j, point[i].y + width / 2 + dy * j);
                if (j < p.nCoils - 1) {
                    ctx.lineTo(point[i].x + dx * (j + 1), point[i].y - width / 2 + dy * j);
                }
                if (j == p.nCoils - 1) {
                    ctx.lineTo(point[i].x + dx * (j + 1), point[i].y - width / 2 + dy * j);
                    ctx.lineTo(point[i + 1].x, point[i + 1].y);
                }
                ctx.stroke();
            }
        }

        if (p.showMass == true) {
            if (i == 0) continue;
            ctx.beginPath();
            ctx.arc(point[i].x, point[i].y, 10 / p.SCALE, 0, 2 * PI);
            ctx.stroke();
            ctx.fill();
        }
    }

    if (p.shakeSpeed != 0 && t < p.timeToDrop) {
        ctx.fillRect(point[0].x - p.width / 4, point[0].y, p.width / 2, p.width / 2);
    }
}

document.getElementById("condition").innerHTML = p.condition;
document.getElementById("strength").innerHTML = p.k.toFixed(2);
var t = 0;
setInterval(() => {
    if (p.timeToDrop != Infinity) {
        document.getElementById("time").innerHTML = "t - " + (Math.floor((p.timeToDrop - t) / 100) + 1);
        if (p.timeToDrop - t < 0) {
            document.getElementById("time").innerHTML = "DROP";
        }
    }
    t += 1;
    if (t % TIME_SLOWINATOR == 0) {
        resetCanvas();
        physics();
        movePoint();
        drawSprings();
    }
}, 1000 * dt);
