var randomNum = (Math.random() * RANDOM_NUM_RANGE - RANDOM_NUM_RANGE / 2).toFixed(2);

var expression = [
    randomNum,
    "x",
    "Math.cos(x)",
    "Math.sin(x)",
    "Math.tan(x)",
    "Math.exp(x)",
    "Math.log(x)",
    "Math.abs(x)"
];
var operator = [
    " + ",
    " * ",
    " ** ",
    " / ",
    " - "
];
var fString = "return ";

function random(mathThing) {
    return mathThing[Math.floor(Math.random() * mathThing.length)]
}


fString += random(expression);
fString += random(operator);
fString += random(expression);
fString += ";";
document.getElementById("function").innerHTML = "f(x) = " + fString.slice(7).replace(/Math./g, "");

var f = new Function(
    "x",
    //"return 0"
    fString
    //"return Math.abs(x / 2 + 0.01) ** Math.cos(x / 2)"
    //"return 1.1 * ((x / 3)**2) ** Math.cos(x / 2)"
    //"return 1.1 * ((x / 2 - 0.03)**2) ** (Math.cos(x / 2) / 1.5)"
)

if (f(Math.random()) == 0 || f(Math.random()) == 1) {
    f = function (x) {
        return 1.1 * ((x / 2 - 0.03) ** 2) ** (Math.cos(x / 2) / 1.5)
    }
    document.getElementById("function").innerHTML = "f(x) = Mustachio de Murray";
}
