#include "Clothoid.hh"
#include "ClothoidAsyPlot.hh"

AsyPlot::AsyPlot(string filename, bool showAxes):
    filename(filename),
    showAxes(showAxes)
{
    if (!openFile()) {
        std::cerr << "Failed to open file " << filename << std::endl;
    }
    initFile();
}

AsyPlot::~AsyPlot() {
    if (showAxes) {
        displayAxes();
    }
    if (closeFile()) {
        compileFile();
    }
}

void AsyPlot::compileFile() {
    string cmdComp = "asy -f pdf "+ filename;
    system(cmdComp.c_str());
    int pos = filename.find(".asy");
    string pdfFile = filename.substr(0,pos) + ".pdf";
    std::cout << pdfFile << std::endl;
    string cmdOpen = "(okular " + pdfFile + " &> /dev/null )&";
    system(cmdOpen.c_str());
}

void AsyPlot::initFile() {
    file << "// File generated automatically from C++ \n\n\n";
    file << "import graph; \n";
    file << "include \"clothoidLib.asylib\"; \n";
    file << "size(14cm,7cm); \n";
    file << "\n\n\n";
}

void AsyPlot::drawClothoid(const Clothoid::ClothoidCurve &c, string penna, double offset)
{
	if (offset == 0.) {
		file << "path pclot = clothoidPoints((";
		file << c.getX0() << ",";
		file << c.getY0() << "),";
		file << c.getTheta0() << ",";
		file << c.getKappa() << ",";
		file << c.getKappa_D() << ",";
		file << c.getL() << ",";
		file << "100,0);\n";
		file << "pen penna = " << penna << ";\n";
		file << "draw(pclot, penna);\n\n";
	}
	else {
		file << "path pclot = clothoidOffset((";
		file << c.getX0() << ",";
		file << c.getY0() << "),";
		file << c.getTheta0() << ",";
		file << c.getKappa() << ",";
		file << c.getKappa_D() << ",";
		file << c.getL() << ",";
		file << "100,";
		file << offset;
		file << "); \n";
		file << "pen penna = " << penna << ";\n";
		file << "draw(pclot, penna);\n\n";
	}
}


void AsyPlot::dot(double x, double y, string penna)
{
    file << "dot((" << x << "," << y << ")," << penna << ");\n\n";
}

void AsyPlot::triangle(const Triangle2D &t, string penna)
{
    file << "draw((" << t.x1() << "," << t.y1() << ") -- " <<
            "(" << t.x2() << "," << t.y2() << ") -- " <<
            "(" << t.x3() << "," << t.y3() << ") -- cycle, " << penna << ");\n\n";
}

void AsyPlot::drawRect(double x0, double y0, double x1, double y1, double x2, double y2, double x3, double y3, string penna) {
	file << "fill((" << x0 << "," << y0 << ") -- " <<
		"(" << x1 << "," << y1 << ") -- " <<
		"(" << x2 << "," << y2 << ") -- " <<
		"(" << x3 << "," << y3 << ") -- cycle, " << penna << ");\n\n";
}

bool AsyPlot::openFile()
{
    file.open(filename);
    return file.is_open();
}

bool AsyPlot::closeFile()
{
    file.close();
    return !file.is_open();
}

void AsyPlot::displayAxes()
{
    file << "xaxis(\"$x$\", black+fontsize(7pt),xmin=0,xmax=5,Ticks(Step=2,step=1,NoZero,Size=.8mm, size=.4mm));\n";
    file << "yaxis(\"$y$\", black+fontsize(7pt),ymin=0,ymax=5,Ticks(Step=2,step=1,NoZero,Size=.8mm, size=.4mm));\n";
}

void AsyPlot::displayAxes(string labX, string labY, double xmin, double xmax, double ymin, double ymax)
{
	file << "xaxis(\"" <<labX<< "\", black+fontsize(7pt),xmin="<<xmin<<",xmax="<<xmax<< ",Ticks(Step=2,step=1,NoZero,Size=.8mm, size=.4mm));\n";
	file << "yaxis(\"" << labY << "\", black+fontsize(7pt),ymin=" << ymin<<",ymax="<<ymax<<",Ticks(Step=2,step=1,NoZero,Size=.8mm, size=.4mm));\n";
}

void AsyPlot::drawLine(double x0, double y0, double x1, double y1, string penna)
{
    file << "draw((" << x0 << "," << y0 << ") -- " <<
         "(" << x1 << "," << y1 << "), " << penna << ");\n\n";
}




void AsyPlot::label(string text, double x, double y, string placement, string penna)
{
	file << "label(\"" << text << "\", (" << x << ", " << y << "), ";
	file << (!placement.empty() ? placement + ", " : "");
	file << penna << ");\n\n";
}


