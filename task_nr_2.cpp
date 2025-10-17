//
// Created by alexf on 15.10.2025.
//
// Graham scan algorithm source: https://www.geeksforgeeks.org/dsa/convex-hull-using-graham-scan
// Closest pair of points algorithm source: https://cp-algorithms.com/geometry/nearest_points.html
//

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
using namespace std;

struct Point {
    float x;
    float y;
};

struct LineFunction {
    float a = 0;
    float b = -1;
    float c = 0;

    float xVal = 0;

    bool isVertical = false;
};

struct Solution3 {
    Point a;
    Point b;
    float distance = FLT_MAX;
};

void printState(vector<Point> &points, int n) {
    for (int i = 0; i < n; i++) {
        cout << i << ". " << to_string(points[i].x) + " " + to_string(points[i].y) << endl;
    }
    cout << endl;
}

// this function loads all the input data
void loadData(vector<Point> &points, int &n) {
    
    cin >> n;

    float x, y;

    for (int i = 0; i < n; i++) {
        cin >> x >> y;
        points.push_back(Point{x, y});
    }

    /*
    n = 10;

    points.push_back({-4, 2});
    points.push_back({2.5, 1.5});
    points.push_back({4, 3});
    points.push_back({0, -2});
    points.push_back({-2, 4});
    points.push_back({4, -1.5});
    points.push_back({2, 2});
    points.push_back({-2, 1});
    points.push_back({-1, 2.5});
    points.push_back({1, 0});
    */
}

// calculates distance using Pythagorean theorem
float distance(Point &a, Point &b) {
    return sqrt ((b.x - a.x) * (b.x - a.x) + (b.y - a.y) * (b.y - a.y));
}

/* 1st part */

// this function returns a value based on if these 3 points make a clockwise turn,
// counterclockwise or are they collinear
// 1 - clockwise
// 0 - collinear
// -1 - counterclockwise
int orientation(Point a, Point b, Point c) {
    //cout << endl << "comparing [" << b.x << ", " << b.y << "] and [" << c.x << ", " << c.y << "]" << endl;
    float crossProduct = ((b.y - a.y) * (c.x - b.x)) - ((b.x - a.x) * (c.y - b.y));

    //cout << "a = [" << a.x << ", " << a.y << "]; b = [" << b.x << ", " << b.y << "]; c = [" << c.x << ", " << c.y << "]\n";

    //cout << crossProduct << endl;

    if (crossProduct == 0) {
        //they are collinear
        return 0;
    }else if (crossProduct < 0) {
        // they are counterclockwise
        return -1;
    }else if (crossProduct > 0) {
        // they are clockwise
        return 1;
    }else {
        cout << "ERROR in orientation" << endl;
        return -2;
    }
}

// this function sorts the points by their polar angle relative to the starting point
void polarAngleSort(vector<Point> &points, int &n) {
    //printState(points, n);

    for (int i = 1; i < n-1;) {

        //cout << "i = " << i << endl;
        //cout << "comparing [" << points[i].x << ", " << points[i].y << "] and [" << points[i+1].x << ", " << points[i+1].y << "]" << endl;

        int orientationVal = orientation(points[0], points[i], points[i+1]);

        //cout << "orientationVal: " << orientationVal << endl;

        if (orientationVal == 0) {
            if (distance(points[0], points[i]) > distance(points[0], points[i+1])) {
                swap(points[i], points[i+1]);
                //cout << "swap" << endl;
            }
            i++;
        }else if (orientationVal == 1) {
            swap(points[i], points[i+1]);
            //cout << "swap" << endl;
            if (i > 0) {
                i--;
            }else {
                i++;
            }
        }else if (orientationVal == -1) {
            //cout << "correct" << endl;
            i++;
        }

        //printState(points, n);

        //cout << endl;

    }

    //cout << "sorted" << endl;
}

// this function finishes the Graham scan algorithm and finds the final list of points
void finishIt(vector<Point> &convexHull, vector<Point> &points, int n) {
    //printState(convexHull, convexHull.size());
    for (int i = 1; i < n; i++) {
        //cout << i << endl;
        //cout << orientation(convexHull[convexHull.size() - 2], convexHull.back(), points[i]) << endl;
        //cout << "[" << points[i].x << "," << points[i].y << "]" << endl;
        //printState(convexHull, convexHull.size());

        while (convexHull.size() > 1 && orientation(convexHull[convexHull.size() - 2], convexHull.back(), points[i]) > 0) {
            convexHull.pop_back();
        }

        convexHull.push_back(points[i]);
    }
}

/*
 * an algorithm for calculating the convex hull of given set of points using Graham scan
 * coded based on description from https://www.geeksforgeeks.org/dsa/convex-hull-using-graham-scan
*/
vector<Point> convexHull(vector<Point> &points, int &n) {
    int minYIndex = 0;
    float minY = points[0].y;
    vector<Point> convexHull = vector<Point>();

    for (int i = 1; i < n; i++) {
        if (points[i].y < minY || (points[i].y == minY && points[i].x < points[minYIndex].x)) {
            minYIndex = i;
            minY = points[i].y;
        }
    }

    convexHull.push_back(points[minYIndex]);

    points[minYIndex] = points[0];
    points[0] = convexHull[0];

    polarAngleSort(points, n);

    finishIt(convexHull, points, n);

    //printState(convexHull, convexHull.size());

    return convexHull;
}

/* 2nd part */

/*
 * this function calculates a, b and c coefficients of a linear function that passes through both points
 */
LineFunction calcLineFunction(Point a, Point b) {
    LineFunction result = LineFunction();

    //cout << "[" << a.x << ", " << a.y << "], [" << b.x << ", " << b.y << "]" << endl;

    if (b.x - a.x != 0) {
        result.a = (b.y - a.y) / (b.x - a.x);
        result.c = a.y - (result.a * a.x);
        //cout << "y = " << result.a << "x + " << result.b << endl;
    }else {
        result.isVertical = true;
        result.xVal = a.x;
        //cout << "it's vertical" << endl;
    }

    return result;
}

/*
 * this function finds the distance between a line and a point
 */
float findDistance(LineFunction f, Point a) {
    //cout << "y = " << f.a << "x + " << f.b << endl;
    //cout << f.isVertical << " " << f.xVal << endl;
    //cout << "[" << a.x << ", " << a.y << "]" << endl;
    float result;

    if (!f.isVertical) {
        result = abs(f.a * a.x + f.b * a.y + f.c);
        result = result / sqrt(f.a * f.a + f.b * f.b);
    }else {
        result = abs(f.xVal - a.x);
    }

    //cout << result << endl << endl;

    if (result == 0) {
        return -1;
    }
    return result;
}

/*
 * this function iterates through all triplets of points on the convex hull and finds
 * the two closest parallel lines that fit all the points between them
 */
float minDistance(vector<Point> &convexHull) {
    float result = 0;
    vector<float> distances = vector<float>(convexHull.size(), 0);
    LineFunction f;

    // this block iterates through all the edges except one
    for (int i = 0; i < convexHull.size() - 1; i++) {
        f = calcLineFunction(convexHull[i], convexHull[i+1]);
        for (int j = 0; j < convexHull.size(); j++) {
            if (j == i) {
                j += 2;
            }
            float temp = findDistance(f, convexHull[j]);
            if (temp > distances[i]) {
                distances[i] = temp;
            }
        }
    }

    // this block checks the last edge
    f = calcLineFunction(convexHull.back(), convexHull[0]);
    for (int i = 1; i < convexHull.size() - 1; i++) {
        float temp = findDistance(f, convexHull[i]);
        if (temp > distances[distances.size() - 1]) {
            distances[distances.size() - 1] = temp;
        }
    }


    // this block finds the smallest value
    //cout << endl;
    result = distances[0];

    for (float distance : distances) {
        //cout << distance << endl;
        if (distance < result && distance > 0) {
            result = distance;
        }
    }

    return result;
}

/* 3rd part */

// merge sort for Points
void myMerge(vector<Point> &points, int left, int mid, int right) {
    int numLeft = mid - left + 1;
    int numRight = right - mid;

    vector<Point> Left(numLeft);
    vector<Point> Right(numRight);

    for (int i = 0; i < numLeft; i++) {
        Left[i] = points[left + i];
    }
    for (int i = 0; i < numRight; i++) {
        Right[i] = points[mid + 1 + i];
    }

    int i = 0;
    int j = 0;
    int k = left;

    while (i < numLeft && j < numRight) {
        if (Left[i].x < Right[j].x || (Left[i].x == Right[j].x && Left[i].y < Right[j].y)) {
            points[k] = Left[i];
            i++;
        }else {
            points[k] = Right[j];
            j++;
        }
        k++;
    }

    while (i < numLeft) {
        points[k] = Left[i];
        i++;
        k++;
    }

    while (j < numRight) {
        points[k] = Right[j];
        j++;
        k++;
    }
}

void mergeSort(vector<Point> &points, int left, int right) {
    if (left >= right) {
        return;
    }

    int mid = left + (right - left) / 2;
    mergeSort(points, left, mid);
    mergeSort(points, mid + 1, right);
    myMerge(points, left, mid, right);
}

// this function checks if two given points are closer than current closest pair
void compareSolution(Point &a, Point &b, Solution3 &solution) {
    float dis = distance(a, b);
    if (dis < solution.distance) {
        solution.distance = dis;
        solution.a = a;
        solution.b = b;
    }
}

void recursiveSolution(int left, int right, vector<Point> &points, int n, Solution3 &solution) {
    if (right - left <= 3) {
        for (int i = left; i < right; i++) {
            for (int j = i+1; j < right; j++) {
                compareSolution(points[i], points[j], solution);
            }
        }
        for (int i = left; i < right - 1; i++) {
            if (points[i].y > points[i+1].y) {
                swap(points[i], points[i+1]);
                if (i > 0) {
                    i -= 2;
                }
            }
        }
        return;
    }

    int mid = (left + right) / 2;
    recursiveSolution(left, mid, points, n, solution);
    recursiveSolution(mid, right, points, n, solution);

    myMerge(points, left, mid, right);

    int tsz = 0; // thanks to this in the second loop we check
    //only the points that are already close to midpoint in x-axis
    for (int i = left; i < right; i++) {
        if (abs(points[i].x - points[mid].x) < solution.distance) {
            for (int j = tsz - 1; j >= 0 && points[i].y - points[j].y < solution.distance; j--) {
                compareSolution(points[i], points[j], solution);
            }
            tsz++;
        }
    }
}

/*
 * this function finds the closest pair in the given set of points,
 * this function is based on algorithm found
 * here: https://cp-algorithms.com/geometry/nearest_points.html
 */
Solution3 closestPair(vector<Point> &points, int &n) {
    Solution3 result = Solution3();

    result.a = points[0];
    result.b = points[n - 1];

    mergeSort(points, 0, n - 1);

    recursiveSolution(0, n - 1,points, n, result);

    return result;
}

int main () {
    vector<Point> points = vector<Point>();
    int n = 0;

    loadData(points, n);

    vector<Point> solution1;

    if (n >= 3) {
        solution1 = convexHull(points, n);
    }else {
        solution1 = vector<Point>();
        for (int i = 0; i < n; i++) {
            solution1.push_back(points[i]);
        };
    }

    //printState(solution1, solution1.size());

    cout << "Otoczka: ";
    for (int i = 0; i < solution1.size(); i++) {
        if (i > 0) {
            cout << ", ";
        }
        cout << "(" << solution1[i].x << ", " << solution1[i].y << ")";
    }
    cout << endl;

    float solution2;

    if (n >= 3) {
        solution2 = minDistance(solution1);
    }else {
        solution2 = -1;
    }

    cout << "Proste: d=" << solution2 << endl;

    Solution3 solution3 = closestPair(points, n);

    cout << "Najblizsze Punkty: [(" << solution3.a.x << ", " << solution3.a.y << "), (" << solution3.b.x << ", " << solution3.b.y << ")] d=" << solution3.distance << endl;

    return 0;
}
