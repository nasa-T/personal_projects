var parameters = [
    //0
    {
        nJoints: 50,
        shakeSpeed: 400 * Math.random() + 100,
        timeToDrop: Infinity,
        springLength: 1,
        nCoils: 1,
        k: Math.random() * 300 + 700,
        friction: 0.99,
        oneShake: true,
        start: 900,
        SCALE: canvas.height / 2000,
        width: 200,
        attach: false,
        sideways: false,
        condition: "Slinky; one shake; no drop; vertical friction"
    },
    //1
    {
        nJoints: 50,
        shakeSpeed: 200 + 300 * Math.random(),
        timeToDrop: 1000,
        springLength: 1,
        nCoils: 1,
        k: Math.random() * 300 + 700,
        friction: 0.99,
        oneShake: false,
        start: 900,
        SCALE: canvas.height / 2000,
        width: 200,
        attach: false,
        sideways: false,
        condition: "Slinky; continuous shake, 10 second drop; vertical friction"
    },
    //2
    {
        nJoints: 50,
        shakeSpeed: 0,
        timeToDrop: Infinity,
        springLength: 1,
        nCoils: 1,
        k: Math.random() * 300 + 700,
        friction: 0.99,
        oneShake: false,
        start: 1000,
        SCALE: canvas.height / 2000,
        width: 200,
        attach: false,
        sideways: false,
        condition: "Slinky; 10 second drop; vertical friction"
    },
    //3
    {
        nJoints: 3,
        shakeSpeed: 0,
        timeToDrop: Infinity,
        springLength: 100,
        nCoils: 20,
        k: Math.random() * 900 + 100,
        friction: 1,
        oneShake: false,
        trace: true,
        notTop: {
            vx: 100,
        },
        start: 250,
        SCALE: canvas.height / 500,
        width: 50,
        attach: false,
        sideways: false,
        condition: "Two springs; bottom spring has rightward velocity; no friction"
    },
    //4
    {
        nJoints: 3,
        shakeSpeed: 0,
        timeToDrop: Infinity,
        springLength: 100,
        nCoils: 20,
        k: 100,
        friction: 1,
        oneShake: false,
        trace: true,
        notTop: {
            x: 100,
        },
        start: 250,
        SCALE: canvas.height / 500,
        width: 50,
        attach: false,
        sideways: false,
        condition: "Two springs; bottom spring has right side offset; no friction"
    },
    //5
    {
        nJoints: 2,
        shakeSpeed: 0,
        timeToDrop: Infinity,
        springLength: 50,
        nCoils: 25,
        k: Math.random() * 500 + 500,
        friction: 1,
        oneShake: false,
        trace: true,
        second: {
            vx: 100,
        },
        start: 50,
        SCALE: canvas.height / 100,
        width: 10,
        attach: false,
        sideways: false,
        condition: "One spring; bottom point has rightward velocity; no friction"
    },
    //6
    {
        nJoints: 2,
        shakeSpeed: 0,
        timeToDrop: Infinity,
        springLength: 50,
        nCoils: 25,
        k: 100,
        friction: 1,
        oneShake: false,
        trace: true,
        notTop: {
            x: 10,
        },
        start: 50,
        SCALE: canvas.height / 100,
        width: 10,
        attach: false,
        sideways: false,
        condition: "One spring; bottom end has right side offset; no friction"
    },
    //7
    {
        nJoints: 3,
        shakeSpeed: 0,
        timeToDrop: Infinity,
        springLength: 100,
        nCoils: 20,
        k: 100,
        friction: 1,
        oneShake: false,
        trace: true,
        second: {
            x: 100,
        },
        start: 250,
        SCALE: canvas.height / 500,
        width: 50,
        attach: true,
        sideways: false,
        condition: "Two top/bottom attached springs; middle point has right side offset; no friction"
    },
    //8
    {
        nJoints: 4,
        shakeSpeed: 0,
        timeToDrop: Infinity,
        springLength: 100,
        nCoils: 15,
        k: 10,
        friction: 1,
        oneShake: false,
        trace: false,
        start: 250,
        SCALE: canvas.height / 500,
        width: 50,
        attach: false,
        second: {
            vx: 100,
        },
        third: {
            vx: -100,
        },
        sideways: true,
        showMass: true,
        condition: "Three horizontal, left/right attached springs; middle masses have equal/opposite velocities; no friction"
    },
    //9
    {
        nJoints: 50,
        shakeSpeed: 200 + 300 * Math.random(),
        timeToDrop: Infinity,
        springLength: 1,
        nCoils: 1,
        k: Math.random() * 300 + 700,
        friction: 0.99,
        oneShake: false,
        start: 900,
        SCALE: canvas.height / 2000,
        width: 200,
        attach: false,
        sideways: false,
        condition: "Slinky; continuous shake; no drop; vertical friction"
    },
    //10
    {
        nJoints: 3,
        shakeSpeed: 0,
        timeToDrop: Infinity,
        springLength: 100,
        nCoils: 20,
        k: Math.random() * 900 + 100,
        friction: 1,
        oneShake: false,
        trace: true,
        second: {
            vx: 100,
        },
        third: {
            vx: -100
        },
        start: 250,
        SCALE: canvas.height / 500,
        width: 50,
        attach: false,
        sideways: false,
        condition: "Two springs; bottom and middle points have equal/opposite velocities; no friction"
    },
    //11
    {
        nJoints: 50,
        shakeSpeed: 200 + 300 * Math.random(),
        timeToDrop: 1000,
        springLength: 1,
        nCoils: 1,
        k: Math.random() * 300 + 700,
        friction: 0.99,
        oneShake: true,
        start: 900,
        SCALE: canvas.height / 2000,
        width: 200,
        attach: false,
        sideways: false,
        condition: "Slinky; one shake, 10 second drop; vertical friction"
    },
    //12
    {
        nJoints: 4,
        shakeSpeed: 0,
        timeToDrop: Infinity,
        springLength: 100,
        nCoils: 20,
        k: 100,
        friction: 1,
        oneShake: false,
        trace: true,
        second: {
            x: Math.random() * 200 - 100,
        },
        third: {
            x: Math.random() * 200 - 100,
        },
        start: 250,
        SCALE: canvas.height / 500,
        width: 50,
        attach: true,
        sideways: false,
        condition: "Three top/bottom attached springs; middle points have two random offsets; no friction"
    },

];