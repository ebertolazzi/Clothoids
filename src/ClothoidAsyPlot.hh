
#ifndef ASYPLOT_H
#define ASYPLOT_H

#include <string>
#include <fstream>
#include <iostream>
#include <cmath>
#include "Clothoid.hh"

using std::string;
using std::ofstream;
using Clothoid::Triangle2D;

class AsyPlot
{
public:
    AsyPlot(string filename, bool showAxes);
    ~AsyPlot();

    void drawClothoid(const Clothoid::ClothoidCurve& c, std::string penna="black", double offset = 0.);
    
    void dot(double x, double y, std::string penna="black");
    void triangle(const Triangle2D& t, std::string penna="black");

	void drawRect(double x0, double y0, double x1, double y1, double x2, double y2, double x3, double y3, string penna="black");

    void drawLine(double x0, double y0, double x1, double y1, std::string penna="black");
	void label(string text, double x, double y, string placement = "", string penna = "black");
	void displayAxes(string labX, string labY, double xmin, double xmax, double ymin, double ymax);

private:
    string filename;
    ofstream file;
    bool showAxes;

    bool openFile();
    bool closeFile();
    void initFile();
    void displayAxes();
    void compileFile();
};


#endif 
