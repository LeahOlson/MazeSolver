/*
 * This is an enumerated type to hold the types of cells required to represent a maze.
 */
MazeCellTypes = {
    WALL: '#',
    PASSAGEWAY: ' ',
    SOLUTION: '@',
    FRONTIER: 'F',
    VISITED: 'V',
}

class MazeCell {
    constructor(row, col, type) {
        this.row = row;
        this.col = col;
        this.type = type;
        // only used by Dijkstra's and A*
        this.priority = Infinity;
    }

    // returns a statement of the cell, its type, and location
    // it is similar to adding a string version of an enum in java
    projection() {
        let projection = '';
        if (this.type === MazeCellTypes.WALL) {
            projection = "WALL cell at "
        }
        // a cell doesn't stop being a passageway when it becomes part of a solution.
        else if (this.type === MazeCellTypes.PASSAGEWAY || this.type === MazeCellTypes.SOLUTION) {
            projection = "PASSAGEWAY cell at "
        }
        projection += '[' + this.row + ',' + this.col + ']';
        return projection;
    }
}

class Maze {
    /*
     *this function constructs the maze object, and accepts an argument,
     * plainTextMaze, which is an nxn string representation of a maze
     * with each cell described by single-character MazeCellTypes.
     */


    // creates the maze array
    // each cell of the maze array is an object 2d array of chars
    constructor(plainTextMaze) {
        // split the string into rows
        this.maze = plainTextMaze.split('\n')

        for (let i = 0; i < this.maze.length; i++) {
            // store each row as a char array
            this.maze[i] = this.maze[i].split('');

            for (let j = 0; j < this.maze[i].length; j++) {
                const type = this.maze[i][j];
                this.maze[i][j] = new MazeCell(i, j, type);
            }
        }

        // hardcoded to have the same start and end
        this.start = this.maze[1][0];
        this.destination = this.maze[this.maze.length - 2][this.maze[0].length - 1];
    }

    /*
     * this function determines whether the argument cell meets the destination criteria
     * returns true if the cell is in the same row and column as the destination cell and false otherwise
     */
    destinationPredicate(cell) {
        return this.destination.row === cell.row && this.destination.col === cell.col;
    }

    /*
     * this function returns all of the neighbors of the argument cell (it does not
     * check whether those neighbors have been visited)
     * Returns direction possibilities in counter clockwise order, change to adding the directions
     * in clockwise order because the end is in the bottom right of the maze.
     */
    getNeighbors(cell) {
        const neighbors = [];

        // test top of cell, now right
        if (cell.col + 1 < this.maze[cell.row].length &&
            this.maze[cell.row][cell.col + 1].type === MazeCellTypes.PASSAGEWAY) {
            neighbors.push(this.maze[cell.row][cell.col + 1]);
        }

        // test left cell, now bottom
        if (cell.row + 1 < this.maze.length &&
            this.maze[cell.row + 1][cell.col].type === MazeCellTypes.PASSAGEWAY) {
            neighbors.push(this.maze[cell.row + 1][cell.col]);
        }

        // test bottom of cell, now left
        if (cell.col - 1 >= 0 &&
            this.maze[cell.row][cell.col - 1].type === MazeCellTypes.PASSAGEWAY) {
            neighbors.push(this.maze[cell.row][cell.col - 1])
        }

        // test right cell, now above
        if (cell.row - 1 >= 0 &&
            this.maze[cell.row - 1][cell.col].type === MazeCellTypes.PASSAGEWAY) {
            neighbors.push(this.maze[cell.row - 1][cell.col]);
        }

        return neighbors;
    }


    /*
     * this function uses a breadth first search to solve the maze. When the solution
     * is found, the function modifies the maze by marking each cell of the solution
     * with the type MazeCellTypes.SOLUTION.
     */
    solveMazeBFS() {
        // create the queue to hold the cells we have visited but need
        let current;
        // to return to explore (we will treat the array like a queue)
        const frontier = [];
        frontier.push(this.start);

        // create a set to hold the cells we have visited and add the
        // first element
        const visited = new Set();
        visited.add(this.start.projection())

        // create a map to hold cells to parents, set first element's
        // parents as false (is source cell). Generally, the parents
        // map will have projection values as keys and objects as values.
        // The parents map holds cells to parents (cells visited prior)
        // the keys are projection cells, the object is the normal cells
        // the variable is necessary to keep track of the cells visited
        const parents = new Array();
        parents[this.start.projection()] = false;

        // search and continue searching  while there are still items in the queue
        while (frontier.length >= 1) {

            // get the next element in the queue
            current = frontier.shift();

            // mark the next element as visited
            current.type = MazeCellTypes.VISITED;

            // test to see if it meets the destination criteria
            if (this.destinationPredicate(current)) {
                // we've reached the destination! Awesome!
                break;
            }

            // get the neighbors of the current cell (passageways)
            const neighbors = this.getNeighbors(current);

            // one by one, add neighbors to the queue
            for (let i = 0; i < neighbors.length; i++) {

                const neighbor = neighbors[i].projection();

                // see if we've already visited this cell
                if (!visited.has(neighbor)) {
                    // if we haven't,  add it to the visited set
                    visited.add(neighbor);
                    // add current as the neighbor's parent
                    parents[neighbor] = current;
                    // add the neighbor to the queue
                    frontier.push(neighbors[i])
                    // set the neighbor to have a "frontier" type
                    neighbors[i].type = MazeCellTypes.FRONTIER;
                }
            }
        }

        // backtrack through each cell's parent and set path cells to type
        // solution
        while (current) {
            current.type = MazeCellTypes.SOLUTION;
            current = parents[current.projection()];
        }
    }

    /*
     * this function uses a depth first search to solve the maze. When the solution
     * is found, the function modifies the maze by marking each cell of the solution
     * with the type MazeCellTypes.SOLUTION.
     */
    solveMazeDFS() {
        // create the queue to hold the cells we have visited but need
        let current;
        // to return to explore (we will treat the array like a queue)
        const frontier = new Array();
        frontier.push(this.start);

        // create a set to hold the cells we have visited and add the
        // first element
        const visited = new Set();
        visited.add(this.start.projection())

        // create a map to hold cells to parents, set first element's
        // parents as false (is source cell). Generally, the parents
        // map will have projection values as keys and objects as values.
        const parents = new Array();
        parents[this.start.projection()] = false;

        // search and continue searching  while there are still items in the queue
        while (frontier.length >= 1) {

            // get the next element in the stack
            current = frontier.pop();

            // mark the next element as visited
            current.type = MazeCellTypes.VISITED;

            // test to see if it meets the destination criteria
            if (this.destinationPredicate(current)) {
                // we've reached the destination! Awesome!
                break;
            }

            // get the neighbors of the current cell (passageways)
            const neighbors = this.getNeighbors(current);

            // one by one, add neighbors to the queue
            for (let i = 0; i < neighbors.length; i++) {

                const neighbor = neighbors[i].projection();

                // see if we've already visited this cell
                if (!visited.has(neighbor)) {
                    // if we haven't,  add it to the visited set
                    visited.add(neighbor);
                    // add current as the neighbor's parent
                    parents[neighbor] = current;
                    // add the neighbor to the queue
                    frontier.push(neighbors[i])
                    // set the neighbor to have a "frotier" type
                    neighbors[i].type = MazeCellTypes.FRONTIER;
                }
            }
        }

        // backtrack through each cell's parent and set path cells to type
        // solution
        while (current) {
            current.type = MazeCellTypes.SOLUTION;
            current = parents[current.projection()];
        }
    }

    /*
     * this function uses a Dijkstra's algorithm to solve the maze. When the solution
     * is found, the function modifies the maze by marking each cell of the solution
     * with the type MazeCellTypes.SOLUTION.
     */
    solveMazeDijkstra() {
        function metric(cell) {
            return cell.priority;
        }

        // create the queue to hold the cells we have visited but need
        let current;
        // to return to explore (we will treat the array like a queue)
        const frontier = new PriorityQueue(metric);
        this.start.priority = 0;
        frontier.enqueue(this.start);

        // create a set to hold the cells we have visited and add the
        // first element
        const visited = new Set();
        visited.add(this.start.projection())

        // create a map to hold cells to parents, set first element's
        // parents as false (is source cell). Generally, the parents
        // map will have projection values as keys and objects as values.
        // The parents map holds cells to parents (cells visited prior)
        // the keys are projection cells, the object is the normal cells
        // the variable is necessary to keep track of the cells visited
        const parents = new Array();
        parents[this.start.projection()] = false;

        // search and continue searching  while there are still items in the queue
        while (frontier.peek()) {

            // get the next element in the queue
            current = frontier.dequeue();

            // if we haven't,  add it to the visited set
            visited.add(current);

            // mark the next element as visited
            current.type = MazeCellTypes.VISITED;

            // test to see if it meets the destination criteria
            if (this.destinationPredicate(current)) {
                // we've reached the destination! Awesome!
                break;
            }

            // get the neighbors of the current cell (passageways)
            const neighbors = this.getNeighbors(current);

            // one by one, add neighbors to the queue
            for (let i = 0; i < neighbors.length; i++) {

                const neighbor = neighbors[i].projection();
                let priority = current.priority + 1;

                // see if we've already visited this cell
                if (!visited.has(neighbor) && priority < neighbors[i].priority) {

                    neighbors[i].priority = priority;
                    // add current as the neighbor's parent
                    parents[neighbor] = current;
                    // add the neighbor to the queue
                    frontier.enqueue(neighbors[i]);
                    // set the neighbor to have a "frontier" type
                    neighbors[i].type = MazeCellTypes.FRONTIER;
                }
            }
        }

        // backtrack through each cell's parent and set path cells to type
        // solution
        while (current) {
            current.type = MazeCellTypes.SOLUTION;
            current = parents[current.projection()];
        }
    }

    /*
     * this function uses A* to solve the maze. When the solution
     * is found, the function modifies the maze by marking each cell of the solution
     * with the type MazeCellTypes.SOLUTION.
     */
    solveMazeAStar() {
        function metric(cell) {
            return cell.priority;
        }

        function heuristic(cell, destination) {
            // Manhattan distance
            let distance =  Math.abs(destination.row - cell.row) + Math.abs(destination.col - cell.col)
            return distance;
        }

        // create the queue to hold the cells we have visited but need
        let current;
        // to return to explore (we will treat the array like a queue)
        const frontier = new PriorityQueue(metric);
        this.start.priority = 0;
        frontier.enqueue(this.start);

        // create a set to hold the cells we have visited and add the
        // first element
        const visited = new Set();
        visited.add(this.start.projection())

        // create a map to hold cells to parents, set first element's
        // parents as false (is source cell). Generally, the parents
        // map will have projection values as keys and objects as values.
        // The parents map holds cells to parents (cells visited prior)
        // the keys are projection cells, the object is the normal cells
        // the variable is necessary to keep track of the cells visited
        const parents = new Array();
        parents[this.start.projection()] = false;

        // search and continue searching  while there are still items in the queue
        while (frontier.peek()) {

            // get the next element in the queue
            current = frontier.dequeue();

            // mark the next element as visited
            current.type = MazeCellTypes.VISITED;

            // test to see if it meets the destination criteria
            if (this.destinationPredicate(current)) {
                // we've reached the destination! Awesome!
                break;
            }

            // get the neighbors of the current cell (passageways)
            const neighbors = this.getNeighbors(current);

            // one by one, add neighbors to the queue
            for (let i = 0; i < neighbors.length; i++) {

                const neighbor = neighbors[i].projection();

                // see if we've already visited this cell
                if (!visited.has(neighbor)) {

                    neighbors[i].priority = current.priority + heuristic(neighbors[i], this.destination);
                    // if we haven't,  add it to the visited set
                    visited.add(neighbor);
                    // add current as the neighbor's parent
                    parents[neighbor] = current;
                    // add the neighbor to the queue
                    frontier.enqueue(neighbors[i]);
                    // set the neighbor to have a "frotier" type
                    neighbors[i].type = MazeCellTypes.FRONTIER;
                }
            }
        }

        // backtrack through each cell's parent and set path cells to type
        // solution
        while (current) {
            current.type = MazeCellTypes.SOLUTION;
            current = parents[current.projection()];
        }
    }


    /*
     * this function returns the number of cells that are included in the solution path.
     */
    cellCounts() {
        const counter = [];
        counter['solution'] = 0;
        counter['visited'] = 0;
        counter['frontier'] = 0;
        for (let i = 0; i < this.maze.length; i++) {
            for (let j = 0; j < this.maze[i].length; j++) {
                if (this.maze[i][j].type === MazeCellTypes.SOLUTION) {
                    counter['solution']++;
                }
                if (this.maze[i][j].type === MazeCellTypes.SOLUTION ||
                    this.maze[i][j].type == MazeCellTypes.VISITED) {
                    counter['visited']++;
                }
                if (this.maze[i][j].type === MazeCellTypes.FRONTIER) {
                    counter['frontier']++;
                }
            }
        }
        return counter;
    }

}

// DFS visits more cells before finding a solution than BFS. On average, DFS visits
// almost all the cells while BFS visits approximately half.
// In the average case, BFS is better suited for pathfinding.