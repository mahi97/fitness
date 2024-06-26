document.getElementById('plotButton').addEventListener('click', function() {
    const anglesInput = document.getElementById('angles').value.split(',').map(Number);
    const magnitudesInput = document.getElementById('magnitudes').value.split(',').map(Number);

    const offset = 60;
    const scale = 1;

    function computeMagnitude(angle, extendedAngles, extendedMagnitudes) {
        angle = angle % 360;
        let applicableRegions = [];
        for (let i = 0; i < extendedAngles.length; i++) {
            let [startAngle, endAngle] = extendedAngles[i];
            let magnitude = extendedMagnitudes[i];
            startAngle = (startAngle + 360) % 360;
            endAngle = (endAngle + 360) % 360;
            if (startAngle > endAngle) {
                if (angle >= startAngle || angle <= endAngle) {
                    let regionCenter = (startAngle + endAngle + 360) / 2 % 360;
                    let distanceToCenter = Math.min(Math.abs(regionCenter - angle), 360 - Math.abs(regionCenter - angle)) || 0.0001;
                    applicableRegions.push([magnitude, distanceToCenter]);
                }
            } else {
                if (startAngle <= angle && angle <= endAngle) {
                    let regionCenter = (startAngle + endAngle) / 2;
                    let distanceToCenter = Math.abs(regionCenter - angle) || 0.0001;
                    applicableRegions.push([magnitude, distanceToCenter]);
                }
            }
        }

        if (applicableRegions.length === 1) {
            return applicableRegions[0][0];
        }
        if (applicableRegions.length > 0) {
            let totalDist = applicableRegions.reduce((acc, [_, d]) => acc + (offset + 1 - d), 0);
            let weights = applicableRegions.map(([_, d]) => (offset + 1 - d) / totalDist);
            let weightedMagnitude = applicableRegions.reduce((acc, [m, _], i) => acc + (m * weights[i]), 0);
            return weightedMagnitude;
        }
        return 0;
    }

    const extendedAngles = anglesInput.map(ang => [(ang - offset), (ang + offset)]);
    const extendedMagnitudes = magnitudesInput.map(mag => mag * scale);
    const targetMagnitudes = Array.from({length: 360}, (_, i) => computeMagnitude(i, extendedAngles, extendedMagnitudes));

    function fourierSeries(theta, a0, a1, b1, a2, b2, a3, b3) {
        return a0 + a1 * Math.cos(theta) + b1 * Math.sin(theta) + a2 * Math.cos(2 * theta) + b2 * Math.sin(2 * theta) + a3 * Math.cos(3 * theta) + b3 * Math.sin(3 * theta);
    }

    function errorFunctionFourier(params) {
        const theta = Array.from({length: 360}, (_, i) => i * 2 * Math.PI / 360);
        const r = theta.map(t => fourierSeries(t, ...params));
        return r.reduce((acc, val, i) => acc + Math.pow(val - targetMagnitudes[i], 2), 0) / 360;
    }

    let initialGuessFourier = [0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5];
    let resultFourier = minimize(errorFunctionFourier, initialGuessFourier);
    let optimizedParamsFourier = resultFourier.x;

    function actionMagnitude(theta) {
        return fourierSeries(theta, ...optimizedParamsFourier);
    }

    function plotMagnitudes() {
        const theta = Array.from({length: 360}, (_, i) => i);
        const optimizedMagnitudesFourier = theta.map(t => actionMagnitude(t * 2 * Math.PI / 360));

        const data = [
            {
                x: theta,
                y: targetMagnitudes,
                mode: 'lines',
                name: 'Target Magnitudes'
            },
            {
                x: theta,
                y: optimizedMagnitudesFourier,
                mode: 'lines',
                name: 'Optimized Fourier Series',
                line: {color: 'red'}
            }
        ];

        Plotly.newPlot('magnitudePlot', data, {title: 'Magnitudes'});
    }

    function plotPath(path, angles, pointA, pointB, title, plotId) {
        const pathTrace = {
            x: path.map(p => p[0]),
            y: path.map(p => p[1]),
            mode: 'lines+markers',
            name: 'Path from A to B'
        };
        const pointATrace = {
            x: [pointA[0]],
            y: [pointA[1]],
            mode: 'markers',
            name: 'Point A',
            marker: {color: 'red'}
        };
        const pointBTrace = {
            x: [pointB[0]],
            y: [pointB[1]],
            mode: 'markers',
            name: 'Point B',
            marker: {color: 'green'}
        };
        const data = [pathTrace, pointATrace, pointBTrace];

        Plotly.newPlot(plotId, data, {title});
    }

    plotMagnitudes();

    const pointA = [28, 27];
    const pointB = [5, 30];

    const { path: pathGreedy, angles: anglesGreedy } = navigateGreedy(pointA, pointB);
    plotPath(pathGreedy, anglesGreedy, pointA, pointB, 'Greedy Path from Point A to Point B using Fourier Series', 'pathPlot');

    const { path: pathAStar, angles: anglesAStar } = aStarSearch(pointA, pointB);
    if (pathAStar.length === 0) {
        console.log("No path found using A*.");
    } else {
        plotPath(pathAStar, anglesAStar, pointA, pointB, 'A* Path from Point A to Point B using Fourier Series', 'pathPlot');
    }
});

// Error minimization function using Nelder-Mead algorithm (simplified version)
function minimize(func, initialGuess) {
    // Placeholder for the actual minimization process
    // This is just a simplified version for demonstration purposes
    let params = initialGuess;
    for (let i = 0; i < 1000; i++) {
        let grad = numericalGradient(func, params);
        params = params.map((p, idx) => p - 0.01 * grad[idx]);
    }
    return { x: params };
}

function numericalGradient(func, params) {
    const h = 1e-5;
    let grad = [];
    for (let i = 0; i < params.length; i++) {
        let params1 = [...params];
        let params2 = [...params];
        params1[i] += h;
        params2[i] -= h;
        grad.push((func(params1) - func(params2)) / (2 * h));
    }
    return grad;
}

function fourierSeries(theta, a0, a1, b1, a2, b2, a3, b3) {
    return a0 + a1 * Math.cos(theta) + b1 * Math.sin(theta) + a2 * Math.cos(2 * theta) + b2 * Math.sin(2 * theta) + a3 * Math.cos(3 * theta) + b3 * Math.sin(3 * theta);
}

function errorFunctionFourier(params) {
    const theta = Array.from({length: 360}, (_, i) => i * 2 * Math.PI / 360);
    const r = theta.map(t => fourierSeries(t, ...params));
    return r.reduce((acc, val, i) => acc + Math.pow(val - targetMagnitudes[i], 2), 0) / 360;
}

function computeMagnitude(angle, extendedAngles, extendedMagnitudes) {
    const offset = 60;
    angle = angle % 360;
    let applicableRegions = [];
    for (let i = 0; i < extendedAngles.length; i++) {
        let [startAngle, endAngle] = extendedAngles[i];
        let magnitude = extendedMagnitudes[i];
        startAngle = (startAngle + 360) % 360;
        endAngle = (endAngle + 360) % 360;
        if (startAngle > endAngle) {
            if (angle >= startAngle || angle <= endAngle) {
                let regionCenter = (startAngle + endAngle + 360) / 2 % 360;
                let distanceToCenter = Math.min(Math.abs(regionCenter - angle), 360 - Math.abs(regionCenter - angle)) || 0.0001;
                applicableRegions.push([magnitude, distanceToCenter]);
            }
        } else {
            if (startAngle <= angle && angle <= endAngle) {
                let regionCenter = (startAngle + endAngle) / 2;
                let distanceToCenter = Math.abs(regionCenter - angle) || 0.0001;
                applicableRegions.push([magnitude, distanceToCenter]);
            }
        }
    }

    if (applicableRegions.length === 1) {
        return applicableRegions[0][0];
    }
    if (applicableRegions.length > 0) {
        let totalDist = applicableRegions.reduce((acc, [_, d]) => acc + (offset + 1 - d), 0);
        let weights = applicableRegions.map(([_, d]) => (offset + 1 - d) / totalDist);
        let weightedMagnitude = applicableRegions.reduce((acc, [m, _], i) => acc + (m * weights[i]), 0);
        return weightedMagnitude;
    }
    return 0;
}

function navigateGreedy(start, end, tol = 0.5) {
    let path = [start];
    let current_position = start;
    let angles = [];
    while (!isPointInside(end, current_position)) {
        let [next_position, best_angle, best_distance] = nextStep(current_position, end);
        path.push(next_position);
        angles.push(best_angle);
        current_position = next_position;
    }
    let [next_position, best_angle, best_distance] = nextStep(current_position, end);
    path.push(next_position);
    angles.push(best_angle);
    return { path, angles };
}

function nextStep(current_position, target_position) {
    let best_distance = Infinity;
    let best_position = current_position;
    let best_angle = 0;
    for (let angle = 0; angle < 2 * Math.PI; angle += Math.PI / 180) {
        let step_size = actionMagnitude(angle);
        let new_position = [
            current_position[0] + step_size * Math.cos(angle),
            current_position[1] + step_size * Math.sin(angle)
        ];
        let distance = Math.hypot(new_position[0] - target_position[0], new_position[1] - target_position[1]);
        if (distance < best_distance) {
            best_distance = distance;
            best_position = new_position;
            best_angle = angle;
        }
    }
    return [best_position, best_angle, best_distance];
}

function isPointInside(target, point) {
    let [x, y] = point;
    let theta = Math.atan2(y, x);
    let r = Math.hypot(x, y);
    let magnitude = actionMagnitude(theta);
    return r <= magnitude;
}

function aStarSearch(start, end) {
    let open_set = new MinHeap();
    open_set.push({ priority: 0, position: start, cost: 0, angle: null });
    let came_from = {};
    let cost_so_far = {};
    cost_so_far[start] = 0;
    let angles = [];
    while (!open_set.isEmpty()) {
        let { position: current, cost: current_cost, angle: current_angle } = open_set.pop();
        if (isPointInside(end, current)) {
            let path = [];
            while (current) {
                path.push(current);
                current = came_from[current];
            }
            path.reverse();
            return { path, angles };
        }
        let possible_moves = getPossibleMovements(current, end);
        for (let [new_position, angle, step_cost] of possible_moves) {
            let new_cost = current_cost + step_cost;
            if (!cost_so_far[new_position] || new_cost < cost_so_far[new_position]) {
                cost_so_far[new_position] = new_cost;
                let priority = new_cost + Math.hypot(new_position[0] - end[0], new_position[1] - end[1]);
                open_set.push({ priority, position: new_position, cost: new_cost, angle });
                came_from[new_position] = current;
                angles.push(angle);
            }
        }
    }
    return { path: [], angles };
}

function getPossibleMovements(current_position, target_position) {
    let movements = [];
    let current_distance = Math.hypot(current_position[0] - target_position[0], current_position[1] - target_position[1]);
    for (let i = 0; i < 36; i++) {
        let angle = 10 * i * Math.PI / 180;
        let step_size = actionMagnitude(angle);
        let new_position = [
            current_position[0] + step_size * Math.cos(angle),
            current_position[1] + step_size * Math.sin(angle)
        ];
        let new_distance = Math.hypot(new_position[0] - target_position[0], new_position[1] - target_position[1]);
        if (new_distance < current_distance) {
            movements.push([new_position, angle, new_distance]);
        }
    }
    movements.sort((a, b) => a[2] - b[2]);
    return movements.slice(0, 2);
}

function MinHeap() {
    this.heap = [];
}

MinHeap.prototype.push = function (item) {
    this.heap.push(item);
    this.heapifyUp(this.heap.length - 1);
};

MinHeap.prototype.pop = function () {
    if (this.heap.length === 1) {
        return this.heap.pop();
    }
    const item = this.heap[0];
    this.heap[0] = this.heap.pop();
    this.heapifyDown(0);
    return item;
};

MinHeap.prototype.isEmpty = function () {
    return this.heap.length === 0;
};

MinHeap.prototype.heapifyUp = function (index) {
    let parent = Math.floor((index - 1) / 2);
    if (index > 0 && this.heap[index].priority < this.heap[parent].priority) {
        [this.heap[index], this.heap[parent]] = [this.heap[parent], this.heap[index]];
        this.heapifyUp(parent);
    }
};

MinHeap.prototype.heapifyDown = function (index) {
    let left = 2 * index + 1;
    let right = 2 * index + 2;
    let smallest = index;
    if (left < this.heap.length && this.heap[left].priority < this.heap[smallest].priority) {
        smallest = left;
    }
    if (right < this.heap.length && this.heap[right].priority < this.heap[smallest].priority) {
        smallest = right;
    }
    if (smallest !== index) {
        [this.heap[index], this.heap[smallest]] = [this.heap[smallest], this.heap[index]];
        this.heapifyDown(smallest);
    }
}
