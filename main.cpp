#include <GL/glut.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
using namespace std;
const int WINDOW_WIDTH = 800;
const int WINDOW_HEIGHT = 800;
double pi = 3.14159265;
bool Rotate = false;
bool Translate = false;
bool shading = false;
double rot_inc=pi/32.0;
int subdivs=0;
double dx=0.0;
double dy=0.0;
double rot_x=0.0;
double rot_y=0.0;
double K_a = .8;
double L_a = .6;
double K_d = .6;
double L_d = .6;
double K_s = .4;
double L_s = .5;
double shiny = 7;


class vertex
{
	public:
	double x,y,z;
	double odd;
	vertex(double s=0,double t=0,double u=0,double o=0)
	:x(s),y(t),z(u),odd(o)
	{}
	friend class Polygon_container;
	vertex operator=(vertex& lhs)
    {
      x=lhs.x;
      y=lhs.y;
      z=lhs.z;
      return *this;
    }
    void print()
    {
		cout << x << " " << y << " " << z << endl;
	}
};

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++
vertex view_orient(-0.2,0.3,0.4);
vertex light_source(0.0,1.0,0.0);
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//overloaded operators for vertex
bool operator==(vertex& rhs, vertex& lhs)
{
  return rhs.x==lhs.x && rhs.y==lhs.y && rhs.z==lhs.z;
};

bool operator!=(vertex& rhs, vertex& lhs)
{
  return !(rhs.x==lhs.x && rhs.y==lhs.y && rhs.z==lhs.z);
};

vertex operator *=(vertex& rhs, double i)
{
	rhs.x*=i;
	rhs.y*=i;
	rhs.z*=i;
	return rhs;
}

vertex operator * ( double i, vertex& rhs)
{
	return vertex(i*rhs.x,i*rhs.y,i*rhs.z);
}

vertex operator += (vertex& rhs, vertex& lhs)
{
	rhs.x+=lhs.x;
	rhs.y+=lhs.y;
	rhs.z+=lhs.z;
	return rhs;
}
vertex operator +(vertex& rhs, vertex& lhs)
{
	return vertex(rhs.x+lhs.x,rhs.y+lhs.y,rhs.z+lhs.z);
}

class Ray
    {
      public:
      vertex pt;
      vertex dir;
      Ray(vertex v,vertex d)
      :pt(v), dir(d)
      {}

};
/**
 * Class that holds all the vertices of a given face inside the polygon.
 * Face is comprised of 3 vertices as polygon is subdivided into triangles.
 */
class face
{
	public:
	vertex a;
	vertex b;
	vertex c;
	face(vertex a, vertex b, vertex c)
	:a(a),b(b),c(c)
	{}
	bool containsVertex(vertex v)
	{
	  if (v==a or v==b or v==c)
	  {
		  return true;
	  }
	  return false;
	}
	double det_t(Ray r)
	{
	  double b_ax = -1*(b.x-a.x);
	  double b_ay = -1*(b.y-a.y);
	  double b_az = -1*(b.z-a.z);
	  double c_ax = -1*(c.x-a.x);
	  double c_ay = -1*(c.y-a.y);
	  double c_az = -1*(c.z-a.z);
	  double a_rx = a.x-r.pt.x;
	  double a_ry = a.y-r.pt.y;
	  double a_rz = a.z-r.pt.z;
	  double det = b_ax*( (c_ay*a_rz)-(c_az*a_ry) ) - (b_ay*( (c_ax*a_rz)-(c_az*a_rx) ) ) + b_az*( (c_ax*a_ry)-(c_ay*a_rx) );
	  return det;
	}
	double det_c(Ray r)
	{
	  double b_ax = -1*(b.x-a.x);
	  double b_ay = -1*(b.y-a.y);
	  double b_az = -1*(b.z-a.z);
	  double a_rx = a.x-r.pt.x;
	  double a_ry = a.y-r.pt.y;
	  double a_rz = a.z-r.pt.z;
	  double det = b_ax*( (a_ry*r.dir.z)-(a_rz*r.dir.y) ) - (b_ay*( (a_rx*r.dir.z)-(a_rz*r.dir.x) ) ) + b_az*( (a_rx*r.dir.y)-(a_ry*r.dir.x) );
	  return det;
	}
	double det_b(Ray r)
	{
	  double a_rx = a.x-r.pt.x;
	  double a_ry = a.y-r.pt.y;
	  double a_rz = a.z-r.pt.z;
	  double c_ax = -1*(c.x-a.x);
	  double c_ay = -1*(c.y-a.y);
	  double c_az = -1*(c.z-a.z);
	  double det = a_rx*( (c_ay*r.dir.z)-(c_az*r.dir.y) ) - (a_ry* ( (c_ax*r.dir.z)-(c_az*r.dir.x) ) ) + a_rz*( (c_ax*r.dir.y)-(c_ay*r.dir.x) );
	  return det;
	}
	double det_a(Ray r)
	{
	  double b_ax = -1*(b.x-a.x);
	  double b_ay = -1*(b.y-a.y);
	  double b_az = -1*(b.z-a.z);
	  double c_ax = -1*(c.x-a.x);
	  double c_ay = -1*(c.y-a.y);
	  double c_az = -1*(c.z-a.z);
	  double det = b_ax*( (c_ay*r.dir.z)-(c_az*r.dir.y) ) - (b_ay*( (c_ax*r.dir.z)-(c_az*r.dir.x) ) ) + b_az*( (c_ax*r.dir.y)-(c_ay*r.dir.x) );
	  return det;
	}
	double get_t(Ray r)
	{
		if (det_a(r)==0)
		{
			return 0.0;
		}
		return det_t(r)/det_a(r);
	}
	double get_B(Ray r)
	{
		if (det_a(r)==0)
		{
			return 0.0;
		}
		return det_b(r)/det_a(r);
	}
	double get_G(Ray r)
	{
		if (det_a(r)==0)
		{
			return 0.0;
		}
		return det_c(r)/det_a(r);
	}
};

//overloaded operators for face
bool operator==(face& rhs, face& lhs)
{
  if (rhs.a != lhs.a){return false;}
  else if (rhs.a != lhs.b){return false;}
  else if (rhs.a != lhs.c){return false;}
  else if (rhs.b != lhs.a){return false;}
  else if (rhs.b != lhs.b){return false;}
  else if (rhs.b != lhs.c){return false;}
  else if (rhs.c != lhs.a){return false;}
  else if (rhs.c != lhs.b){return false;}
  else if (rhs.c != lhs.c){return false;}
  else{return true;}
};


//comparator for map objects used in polygon container
struct compare_v
{
   bool operator() (const vertex& v, const vertex& w)
   {
          if (v.x==w.x)
          {
                  if (v.y == w.y)
                  return v.z<w.z;
                  else
                  return v.y < w.y;
          }
          return v.x < w.x;
   }
};

/**
 * Class that stores and contains all the data used for the Loop subdivision.
 * Container starts with a set of vertices that are then averaged after subdividing.
 * New faces are also generated off of the original faces after Loop subdivision has
 * been applied.
 */
class Polygon_container
{
	public:
	vector<vertex> verts;
	vector<vertex> scaled_verts;
	vector<vertex> trans_scaled_verts;
	vector<face> faces;
	vector <face> new_faces;
	vector<face> trans_scaled_faces;
	vector<face> scaled_new_faces;
	vector<vertex>norms;
	map< int,vector<int> > face_rel;
    map<vertex,vector<vertex>,compare_v> neighbors;
	Polygon_container()
	{}
	bool isNeighbor(vertex home, vertex neighbor)
	{
		if (neighbors[home].size()==0)
		return false;
		for (int i(0); i < neighbors[home].size(); ++i)
	    {
		   if (neighbors[home][i]==neighbor)
		   return true;
		}
		return false;
	}

	//finds a vertex in the verts vector
	int whereIs(vertex& v)
	{
		int counter(0);
		for(; counter < verts.size(); ++counter)
		{
			if (verts[counter]==v)
			{
				return counter;
			}
		}
	}
	int whereIs1(vector<vertex>& t,vertex& v)
	{
		int counter(0);
		for(; counter < t.size(); ++counter)
		{
			if (t[counter]==v)
			{
				return counter;
			}
		}
	}

	//finds center of polygon
	vertex center()
	{
	   double sumx(0.0),sumy(0.0),sumz(0.0);
	   for (int i(0); i < trans_scaled_verts.size(); ++i)
	   {
		  sumx+=trans_scaled_verts[i].x;
		  sumy+=trans_scaled_verts[i].y;
		  sumz+=trans_scaled_verts[i].z;
	   }
	   vertex cent(sumx/double(scaled_verts.size()),sumy/double(scaled_verts.size()),sumz/double(scaled_verts.size()));
	   return cent;
	}

    //checks whether a vertex a is in vector of vertices
	bool isIn(vertex a, vector<vertex>& v)
	{
		for (int i(0); i < v.size(); ++i)
		{
		   if (v[i]==a)
		   {
		     return true;
		   }
		}
		return false;
	}

	//generates the neighbors for the polygon at its initial state
	void generate_neighbors()
	{
	  for (int i(0); i < faces.size(); i++)
	  {
		  for (int j(0); j < 3; ++j)
		  {
               if (j==0)
               {
				   if(!isNeighbor(faces[i].a,faces[i].b))
				   {
				     neighbors[faces[i].a].push_back(faces[i].b);
				     neighbors[faces[i].b].push_back(faces[i].a);
				   }
				   if (!isNeighbor(faces[i].a,faces[i].c))
				   {
				     neighbors[faces[i].a].push_back(faces[i].c);
				     neighbors[faces[i].c].push_back(faces[i].a);
				   }
			   }
			   else if (j==1)
			   {
			     if(!isNeighbor(faces[i].b,faces[i].a))
				 {
				    neighbors[faces[i].b].push_back(faces[i].a);
				    neighbors[faces[i].a].push_back(faces[i].b);
				 }
				 if (!isNeighbor(faces[i].b,faces[i].c))
				 {
				    neighbors[faces[i].b].push_back(faces[i].c);
				    neighbors[faces[i].c].push_back(faces[i].b);
				 }
			   }
			   else
			   {
			     if(!isNeighbor(faces[i].c,faces[i].a))
				 {
				    neighbors[faces[i].c].push_back(faces[i].a);
				    neighbors[faces[i].a].push_back(faces[i].c);
				 }
				 if (!isNeighbor(faces[i].c,faces[i].b))
				 {
				   neighbors[faces[i].c].push_back(faces[i].b);
				   neighbors[faces[i].b].push_back(faces[i].c);
				 }
			   }
		  }
	  }
	}

	//generates the vertices neighbor relationship amonst the other vertices in the
	//polygon
	void generate_neighbor_rel()
	{
	  for (int i(0); i < verts.size();++i)
	  {
		 for(int j(0); j < neighbors[verts[i]].size(); ++j)
		 {
			 int ind = whereIs(neighbors[verts[i]][j]);
			 face_rel[i].push_back(ind);
		 }
	  }
    }

    //populates the scaled faces vector and is called after Loop subdivision is da
    void generate_scaled_faces()
    {
		for (int i(0); i < new_faces.size(); ++i)
		{
			face f = new_faces[i];
			int a = whereIs(new_faces[i].a);
			int b = whereIs(new_faces[i].b);
			int c = whereIs(new_faces[i].c);
			vertex s = scaled_verts[a];
			vertex t = scaled_verts[b];
			vertex u = scaled_verts[c];
			face scale(s,t,u);
			scaled_new_faces.push_back(scale);
		}

	}

	//generates the faces after transformations have been applied
	void generate_trans_scaled_faces()
	{
		trans_scaled_faces.clear();
		for (int i(0); i < scaled_new_faces.size(); ++i)
		{
		  	face f = scaled_new_faces[i];
			int a = whereIs1(scaled_verts,f.a);
			int b = whereIs1(scaled_verts,f.b);
			int c = whereIs1(scaled_verts,f.c);
			vertex s = trans_scaled_verts[a];
			vertex t = trans_scaled_verts[b];
			vertex u = trans_scaled_verts[c];
			face scale(s,t,u);
			trans_scaled_faces.push_back(scale);
		}
	}

	//rotates the polygon
	void rotate()
	{
		vertex cent = center();
		for (int i(0); i < trans_scaled_verts.size(); ++i)
		{
		  trans_scaled_verts[i].x-=cent.x;
		  trans_scaled_verts[i].y-=cent.y;
		  trans_scaled_verts[i].z-=cent.z;
		}
		for (int j(0); j < trans_scaled_verts.size(); ++j)
		{
			double temp_x = trans_scaled_verts[j].x;
			double temp_y=trans_scaled_verts[j].y;
		   	trans_scaled_verts[j].x=temp_x*cos(rot_x)-trans_scaled_verts[j].z*sin(rot_x);
		   	trans_scaled_verts[j].z=trans_scaled_verts[j].z*cos(rot_x)+temp_x*sin(rot_x);
		   	trans_scaled_verts[j].y=temp_y*cos(rot_y)-trans_scaled_verts[j].z*sin(rot_y);
		   	trans_scaled_verts[j].z=trans_scaled_verts[j].z*cos(rot_y)+temp_y*sin(rot_y);
		}
		for (int o(0); o < trans_scaled_verts.size(); ++o)
		{
		  trans_scaled_verts[o].x+=cent.x;
		  trans_scaled_verts[o].y+=cent.y;
		  trans_scaled_verts[o].z+=cent.z;
		}

	}

	//creates all the vertices that are translated from original polygon mesh
	void generate_trans_scaled_verts()
	{
		trans_scaled_verts.clear();
		for (int i(0); i < scaled_verts.size(); ++i)
		{
		   vertex v(scaled_verts[i].x+dx,scaled_verts[i].y+dy,scaled_verts[i].z);
		   trans_scaled_verts.push_back(v);
		}
		rotate();
	}

	//does cross product between 2 "vectors", I used the vertex class to represent
	//a linear algebra vector
	vertex cross_prod(vertex& v1, vertex&v2)
	{
	  	vertex n;
	  	n.x=(v1.y*v2.z)-(v1.z*v2.y);
	  	n.y=(v1.z*v2.x)-(v1.x*v2.z);
	  	n.z=(v1.x*v2.y)-(v1.y*v2.x);
	  	return n;
	}

	//does the dota product between 2 "vectors", I used the vertex class to represent
	//a linear algebra vector
	double dot_prod(vertex& v1, vertex& v2)
	{
		return v1.x*v2.x+v1.y*v2.y+v1.z*v2.z;
	}
	//normalizes a normal to make a unit "vector", I used the vertex class to represent
	//a linear algebra vector
	void normalize(vertex& norm)
	{
	   double div = sqrt(norm.x*norm.x+norm.y*norm.y+norm.z*norm.z);
	   norm.x/=div;
	   norm.y/=div;
	   norm.z/=div;
	}

	//calculates the normalized norms of each face with respect to the polygon's
	//center and stores them in a vector
    void calc_norms()
    {
		norms.clear();
	   	for (int i(0); i < trans_scaled_faces.size();++i)
	   	{
			vertex c = center();
			vertex v1(trans_scaled_faces[i].b.x-trans_scaled_faces[i].a.x,trans_scaled_faces[i].b.y-trans_scaled_faces[i].a.y,trans_scaled_faces[i].b.z-trans_scaled_faces[i].a.z);
			vertex v2(trans_scaled_faces[i].c.x-trans_scaled_faces[i].a.x,trans_scaled_faces[i].c.y-trans_scaled_faces[i].a.y,trans_scaled_faces[i].c.z-trans_scaled_faces[i].a.z);
			vertex norm = cross_prod(v1,v2);
			normalize(norm);
			vertex new_vert(c.x-trans_scaled_faces[i].b.x,c.y-trans_scaled_faces[i].b.y,c.z-trans_scaled_faces[i].b.z);
			if (dot_prod(norm,new_vert)>0.0)
			{
				norm.x*=-1;
				norm.y*=-1;
				norm.z*=-1;
			}
			norms.push_back(norm);
		}
	}
    //Loops subdivision: divides each face into 4 triangular faces and sets up the process for
    //the smoothing/averaging part to happen
	void subdivide()
	{
		for (int i(0); i < faces.size(); ++i)
		{
		   	vertex a = faces[i].a;
		   	vertex b = faces[i].b;
		   	vertex c = faces[i].c;
		   	vertex mid1((a.x+b.x)/2,(a.y+b.y)/2,(a.z+b.z)/2,1);
		   	vertex mid2((a.x+c.x)/2,(a.y+c.y)/2,(a.z+c.z)/2,1);
		   	vertex mid3((c.x+b.x)/2,(c.y+b.y)/2,(c.z+b.z)/2,1);
		   	if (!isIn(mid1,verts))
		   	{
		   	  verts.push_back(mid1);
		   	}
		   	if (!isIn(mid2,verts))
		   	{
		   	  verts.push_back(mid2);
		   	}
		   	if (!isIn(mid3,verts))
		   	{
		   	  verts.push_back(mid3);
		   	}
		   	face face1(a,mid1,mid2);
		   	face face2(b,mid1,mid3);
		   	face face3(c,mid2,mid3);
		   	face face4(mid1,mid2,mid3);
		   	new_faces.push_back(face1);
		   	new_faces.push_back(face2);
		   	new_faces.push_back(face3);
		   	new_faces.push_back(face4);
		   	if (!isIn(mid1,neighbors[a]))
		   	{
				neighbors[a].push_back(mid1);
				neighbors[mid1].push_back(a);
			}
			if (!isIn(mid2,neighbors[a]))
		   	{

				neighbors[a].push_back(mid2);
				neighbors[mid2].push_back(a);
			}
			if (!isIn(mid1,neighbors[b]))
		   	{
				neighbors[b].push_back(mid1);
				neighbors[mid1].push_back(b);
			}
			if (!isIn(mid3,neighbors[b]))
		   	{
				neighbors[b].push_back(mid3);
				neighbors[mid3].push_back(b);
			}
			if (!isIn(mid2,neighbors[c]))
		   	{
				neighbors[c].push_back(mid2);
				neighbors[mid2].push_back(c);
			}
			if (!isIn(mid3,neighbors[c]))
		   	{
				neighbors[c].push_back(mid3);
				neighbors[mid3].push_back(c);
			}

			if (!isIn(mid1,neighbors[mid2]))
			{
				neighbors[mid2].push_back(mid1);
			}
			if (!isIn(mid1,neighbors[mid3]))
			{
				neighbors[mid3].push_back(mid1);
			}
			if (!isIn(mid2,neighbors[mid1]))
			{
				neighbors[mid1].push_back(mid2);
			}
			if (!isIn(mid2,neighbors[mid3]))
			{
				neighbors[mid3].push_back(mid2);
			}
			if (!isIn(mid3,neighbors[mid1]))
			{
				neighbors[mid1].push_back(mid3);
			}
			if (!isIn(mid3,neighbors[mid2]))
			{
				neighbors[mid2].push_back(mid3);
			}

		}
		generate_neighbor_rel();
	}

	//averages the points in the polygon and generates the scaled faces vector
	void average()
	{
	  for (int i(0); i < verts.size(); ++i)
	  {
		  if (verts[i].odd == 1)
		  {
			  verts[i].odd=0;
			  scaled_verts.push_back(verts[i]);
		  }
		  else
		  {

			  double k = neighbors[verts[i]].size();
			  double c= cos(2*M_PI/k);
			  c*=(1.0/4.0);
			  c+=(3.0/8.0);
			  double d = pow(c,2.0);
			  double B = (1.0/k)*((5.0/8.0)-d );
			  vertex scale = (1-(k*B))*verts[i];
			  for (int j(0); j < neighbors[verts[i]].size(); ++j)
			  {
				  vertex add = neighbors[verts[i]][j];
				  add*=B;
				  scale+=add;
			  }
			  scaled_verts.push_back(scale);
		  }
	  }
	  generate_scaled_faces();
    }

   //used for switching between iterations so that you know when you subdivide the
   //difference between an even vertex and an odd vertex
   void all_even()
   {
	 for (int i(0); i > scaled_verts.size(); ++i)
	 {
	    scaled_verts[i].odd=0;
	 }
	 for ( int j(0); j <  scaled_new_faces.size(); ++j)
	 {
	   	 scaled_new_faces[j].a.odd=0;
	   	 scaled_new_faces[j].b.odd=0;
	   	 scaled_new_faces[j].c.odd=0;
	 }
   }

   //prepares the next iteration's polygon container with the verts and faces that it will use to
   //do Loop's subdivision on
   void transfer_info(Polygon_container& p)
   {
	   	all_even();
	    p.verts.clear();
	    p.faces.clear();
	    copy(scaled_verts.begin(),scaled_verts.end(),back_inserter(p.verts));
	    copy(scaled_new_faces.begin(),scaled_new_faces.end(),back_inserter(p.faces));
   }

   //generates the maxes and minimums of the polygon to calculate the shading faster
   void max_min(int&xmin, int&xmax,int&ymin,int&ymax)
   {
	   xmin=int(trans_scaled_verts[0].x);
	   xmax=int(trans_scaled_verts[0].x);
	   ymin=int(trans_scaled_verts[0].y);
	   ymax=int(trans_scaled_verts[0].y);
	   for (int i(1); i < trans_scaled_verts.size(); ++i)
	   {
		   if (trans_scaled_verts[i].x < xmin)
		   {
			  xmin=int(trans_scaled_verts[i].x)-1;
		   }
		   if (trans_scaled_verts[i].x > xmax)
		   {
			  xmax=int(trans_scaled_verts[i].x)+1;
		   }
		   if (trans_scaled_verts[i].y < ymin)
		   {
			  ymin=int(trans_scaled_verts[i].y)-1;
		   }
		   if (trans_scaled_verts[i].y > ymax)
		   {
			  ymax=int(trans_scaled_verts[i].y)+1;
		   }
	   }
   }
};

Polygon_container zero_it;
Polygon_container first_it;
Polygon_container sec_it;
Polygon_container third_it;
// Renders a quad at cell (x, y) with dimensions CELL_LENGTH
void renderPixel(double x, double y, float r = 1.0, float g = 1.0, float b = 1.0)
{
glBegin(GL_POINTS);
glColor3f(r, g, b);
glVertex2f(x,y);
glEnd();
}

void renderLine(double x0, double y0, double x1, double y1, double r,double g,double b)
{
    glBegin(GL_POINTS);
    if (x0==x1)
    {
		if (y0>y1)
		{
		  double t = y0;
		  y0=y1;
		  y1=t;
		}
		for(;y0<=y1;++y0)
		renderPixel(x0,y0,r,g,b);
	}
	else
	{
	  if (x0 > x1)
	  {
		  double temp = x0;
		  x0=x1;
		  x1=temp;
		  double t = y0;
		  y0=y1;
		  y1=t;
	  }
      double i = x0;
      double m = (y1-y0)/(x1-x0);
      double y = y0;
      while(i<=x1)
      {
		 renderPixel(i,round(y),r,g,b);
		 ++i;
		 y+=m;
      }
     }
    glEnd();
}

//gets the info of the mesh from a text file and stores it in the 1st iterations
//polygon container
void get_info(char** argv, Polygon_container& p)
{
	ifstream read;
	read.open(argv[1]);
	int vert,faces;
	read >> vert >> faces;
	for (int i(0); i < vert; ++i)
	{
	  vertex val;
	  read >> val.x >> val.y >> val.z;
	  val*=8.0;
	  val.odd=0;
	  p.verts.push_back(val);
	}
	for (int j(0); j < faces; ++j)
	{
	  int v1,v2,v3;
	  read >> v1 >> v2 >> v3;
	  face fac(p.verts[v1],p.verts[v2],p.verts[v3]);
	  p.faces.push_back(fac);
	}
	p.generate_neighbors();
    read.close();
}






//gets keyboard input to do transformations, phung shading, and shadow rays
void specialKey( unsigned char key, int x, int y )
{

       switch (key)
       {
		 case 'i':
			{
			   if (shading==true)
			   {
				   shading = false;
			   }
			   else
			   {
				   shading=true;
				   Translate=false;
				   Rotate=false;
			   }
			   break;
			}
		 case 'l':
			{
			   if (subdivs==3)
			   {
				  subdivs=0;
				  break;
			   }
			   else
			   {
				 ++subdivs;
			   }
			   break;
			}
		 case 'r':
			{

			    Translate=false;
				Rotate = true;
				break;
			}
		case 't':
			{

				Rotate=false;
				Translate = true;
				break;
			}
		case 'w':
			{
				if (Rotate and !Translate)
				{
					rot_y += rot_inc;
					break;
				}
				else if (Translate and !Rotate)
				{
					dy+=1.0;
				}
				break;
			}
		case 'a':
			{
				if (Rotate and !Translate)
				{
					rot_x-=rot_inc;
					break;
				}
				else if (Translate and !Rotate)
				{
					dx-=1.0;
				}
				break;
			}
		case 's':
			{
				if (Rotate and !Translate)
				{
					rot_y -= rot_inc;
					break;
				}
				else if (Translate and !Rotate)
				{
					dy -= 1.0;
				}
				break;
			}
		case 'd':
			{
				if (Rotate and !Translate)
				{
					rot_x+=rot_inc;
					break;
				}
				else if (Translate and !Rotate)
				{
					dx+=1.0;
				}
				break;
			}
	   }
	   glutPostRedisplay();
	   return;
}

//used to render polygon's wire fram when shading and shows are not applied
void wire_frame(Polygon_container& p)
{
	 	 p.generate_trans_scaled_verts();
	    for (int i(0); i < p.trans_scaled_verts.size(); ++i)
	    {
		  for (int j(0); j < p.face_rel[i].size(); ++j)
		  {
			  int ind = p.face_rel[i][j];
			  renderLine(p.trans_scaled_verts[i].x,p.trans_scaled_verts[i].y,p.trans_scaled_verts[ind].x,p.trans_scaled_verts[ind].y,0.0,1.0,0.0);
		  }
	    }
}

//globalized dot product of 2 "vectors", I used the vertex class to represent
//a linear algebra vector
double d_prod(vertex& v1, vertex& v2)
{
  return v1.x*v2.x+v1.y*v2.y+v1.z*v2.z;
}

//phong shader that uses constants, light source, and the viewer's orientation
//defined at top of file to shade a triangular mesh polygon
void p_shade(Polygon_container&p)
{
     int xmin = 0,xmax = 0,ymin = 0,ymax = 0;
     p.max_min(xmin,xmax,ymin,ymax);
     for (int i(xmin); i < xmax; ++i)
     {
	     for (int j(ymax); j >ymin; --j)
	     {
            vertex v(i,j,0);
            vertex dir(0,0,1);
		    Ray r(v,dir);
		    double min(232323.0);
		    int ind(-1);
		    p.generate_trans_scaled_faces();
		      for (int k(0); k < p.trans_scaled_faces.size();++k)
		      {
			    if (p.trans_scaled_faces[k].get_t(r)>0)
		        {
		          double b = p.trans_scaled_faces[k].get_B(r);
		          double g = p.trans_scaled_faces[k].get_G(r);
		          if (b > 0 and g > 0)
		          {
			        if (b+g<1)
			        {
						   if (p.trans_scaled_faces[k].get_t(r)<min)
						   {
				              min = p.trans_scaled_faces[k].get_t(r);
				              ind = k;
				            }
			        }
			        else
		            {
				      continue;
				    }
		          }
		        }
		        else
		        {
				  continue;
				}
			  }
			  if (min!=232323.0)
			  {
				 double ambient_light=K_a*L_a;
				 double diffuse_light=K_d*L_d*d_prod(light_source,p.norms[ind]);
				 double dot=(-2.0)*d_prod(light_source,p.norms[i]);
				 vertex ref(ref.x=light_source.x-2.0*d_prod(light_source,p.norms[ind])*p.norms[ind].x, ref.y=light_source.y-2.0*d_prod(light_source,p.norms[ind])*p.norms[ind].y, ref.z=light_source.z-2.0*d_prod(light_source,p.norms[ind])*p.norms[ind].z);
				 double specular_light=K_s*pow(d_prod(ref,light_source),shiny);
				 double color=(ambient_light+diffuse_light+specular_light);
				 renderPixel(i,j,specular_light,specular_light,color/K_a);
			  }
	     }
	 }
}

//creates a shadow on a ground plane by ray tracing from the ground plane to the polygon
//if the ray hits the polygon that pixel is rendered black, else it will be rendered green
//to represent the ground
void shadow(Polygon_container& p)
{
	for (int i(0); i < WINDOW_WIDTH; ++i)
	{
		for (int j(300); j > 0;--j)
		{
			vertex v(i,j,0);
            vertex dir;
            dir.x=view_orient.x;
            dir.y=view_orient.y;
            dir.z=view_orient.z;
		    Ray r(v,dir);
		    double min(232323.0);
		    int ind(-1);
		    p.generate_trans_scaled_faces();
			for (int k(0); k < p.trans_scaled_faces.size();++k)
		    {
			  if (p.trans_scaled_faces[k].get_t(r)>0)
		      {
		         double b = p.trans_scaled_faces[k].get_B(r);
		         double g = p.trans_scaled_faces[k].get_G(r);
		         if (b > 0 and g > 0)
		         {
			       if (b+g<1)
			       {
					   min = p.trans_scaled_faces[k].get_t(r);
			           ind = k;
			       }
			       else
		           {
				     continue;
				   }
		        }
		      }
		      else
		      {
				continue;
			  }
			}
			if (ind==-1)
			{
			  renderPixel(i,j,0.0,1.0,0);
			}
			else
			{
				renderPixel(i,j,0.0,0.0,0.0);
			}
		}
	}
}

void GL_render()
{

   glClear(GL_COLOR_BUFFER_BIT);
   if (shading==false)
   {
	  if (subdivs==0)
	  {
		 zero_it.generate_trans_scaled_verts();
	     wire_frame(zero_it);
	   }
	   else if (subdivs==1)
	   {
		  first_it.generate_trans_scaled_verts();
	      wire_frame(first_it);
	   }
	   else if (subdivs==2)
	   {
		  sec_it.generate_trans_scaled_verts();
	      wire_frame(sec_it);
	   }
	   else
	   {
		  third_it.generate_trans_scaled_verts();
	      wire_frame(third_it);
	   }
   }
   else
   {
	  cout << "shading..."<<endl;
	  if (subdivs == 0)
	  {
		  zero_it.generate_trans_scaled_faces();
		  zero_it.calc_norms();
		  shadow(zero_it);
	      p_shade(zero_it);
	  }
	  else if (subdivs == 1)
	  {
		  first_it.generate_trans_scaled_faces();
		  first_it.calc_norms();
		  shadow(first_it);
	      p_shade(first_it);
	  }
	  else if (subdivs == 2)
	  {
		   sec_it.generate_trans_scaled_faces();
		   sec_it.calc_norms();
		   shadow(sec_it);
	       p_shade(sec_it);
	  }
	  else if (subdivs == 3)
	  {
		   third_it.generate_trans_scaled_faces();
		   third_it.calc_norms();
		   shadow(third_it);
	       p_shade(third_it);
	  }
	  cout << "done shading"<<endl;
	}
   glutSwapBuffers();
}
//Initializes OpenGL attributes
void GLInit(int* argc, char** argv)
{
	glutInit(argc, argv);
	glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE);
	glutInitWindowSize(WINDOW_WIDTH, WINDOW_HEIGHT);

	// ...
	// Complete this function
	// ...
	glutCreateWindow("CS 130 - Tyler Wilson");

	// The default view coordinates is (-1.0, -1.0) bottom left & (1.0, 1.0) top right.
	// This is set to the number of pixels
	// in each dimension.
	glMatrixMode(GL_PROJECTION_MATRIX);
	glOrtho(0, WINDOW_WIDTH, 0, WINDOW_HEIGHT, -1, 1);
    glutDisplayFunc(GL_render);
    glutKeyboardFunc(specialKey);
}

int main(int argc, char** argv)
{
	get_info(argv, zero_it);
	copy(zero_it.verts.begin(),zero_it.verts.end(),back_inserter(zero_it.scaled_verts));
	copy(zero_it.faces.begin(),zero_it.faces.end(),back_inserter(zero_it.new_faces));
    copy(zero_it.faces.begin(),zero_it.faces.end(),back_inserter(zero_it.scaled_new_faces));
    zero_it.generate_neighbor_rel();
    zero_it.transfer_info(first_it);
    first_it.subdivide();
    first_it.average();
    first_it.transfer_info(sec_it);
    sec_it.subdivide();
    sec_it.average();
    sec_it.transfer_info(third_it);
    third_it.subdivide();
    third_it.average();
	GLInit(&argc, argv);
	glutMainLoop();

	return 0;
}
