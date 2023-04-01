#include "FPToolkit.c"
#include "M3d_matrix_tools.c"
#include "xwd_tools_03.c"

double eye[3], coi[3], up[3] ;
double world[3];
int idA, idB ;
int textureWidth, textureHeight;

#include "partials.c"
#include "shape_constructors.c"
#include "light_model.c"
#include "rgb_mapping.c"
#include "ray.c"

void gather_textures(){
  char nameA[100], nameB[100] ;
  int widthA,heightA ;

  int d[2],d1[2],w,h,e,g ;
  w = 1 ; h = 1 ;
  printf("enter name of xwd file\n") ;
  scanf("%s",nameA) ;
  idA = init_xwd_map_from_file (nameA) ;// returns -1 on error, 1 if ok
  if (idA == -1) { printf("failure\n") ;  exit(0) ; }
  e = get_xwd_map_dimensions(idA, d) ;
  if (e == -1) { printf("failure\n") ;  exit(0) ; }
  textureWidth = d[0] ; textureHeight = d[1] ;

  printf("enter name of xwd file\n") ;
  scanf("%s",nameB) ;
  idB = init_xwd_map_from_file (nameB) ;// returns -1 on error, 1 if ok
  if (idB == -1) { printf("failure\n") ;  exit(0) ; }
  g = get_xwd_map_dimensions(idB, d1) ;
  if (g == -1) { printf("failure\n") ;  exit(0) ; }
  textureWidth = d1[0] ; textureHeight = d1[1] ;
}


int main ()
{

    gather_textures();

    int window_width, window_height, window_square_size ;
    double Half_window_size ;
    double Half_angle_degrees;
    double Tan_half_angle ;

    window_width = 800 ; window_height = 800 ;
    // size of largest square INside window
    if (window_width < window_height) { window_square_size = window_width ; }
                                 else { window_square_size = window_height ; }
    Half_window_size = 0.5*window_square_size ;
    Half_angle_degrees = 30 ;
    Tan_half_angle = tan(Half_angle_degrees*M_PI/180) ;

    G_init_graphics (window_width, window_height) ;

  	int q;
    double del = M_PI/50, angle = 0, t = 0;
    double v[3] = {sin(angle), 0, cos(angle)};

    eye[0] = 0 ; eye[1] = del ; eye[2] = 0 ;
    world[0] = 0; world[1] = 0; world[2] = 0;

  	while (q != 'q') {

      coi[0] = eye[0] + v[0];
      coi[1] = eye[1] + v[1];
      coi[2] = eye[2] + v[2];

      up[0] = eye[0] + 0;
      up[1] = eye[1] + 1;
      up[2] = eye[2] + 0;

      M3d_view (vm, vi,  eye, coi, up) ; //create view matrix

      G_rgb(0,0,0) ;
      G_clear() ;

      num_objects = 0;
      color[num_objects][0] = 0.2 ; // grey
      color[num_objects][1] = 0.2 ;
      color[num_objects][2] = 0.2 ;
      ref[num_objects] = 0;
    	sphere(-.5, .1, 1, .1);
      num_objects++;
      color[num_objects][0] = 0.0 ; // blue
      color[num_objects][1] = 0.0 ;
      color[num_objects][2] = 1.0 ;
      ref[num_objects] = .2; // reflects nothing
      sphere(1, .1, .1, .1);
      num_objects++;
      color[num_objects][0] = 1 ; // pink
      color[num_objects][1] = 0.0 ;
      color[num_objects][2] = 1 ;
      ref[num_objects] = .2; // high reflection
      plane(0, -.1, 1, 1, 0);
      num_objects++;
      ref[num_objects] = .2; // high reflection
      plane(0, 0, 2, 1, 90);
      num_objects++;
      color[num_objects][0] = 0.0 ; // blue
      color[num_objects][1] = 0.0 ;
      color[num_objects][2] = 1.0 ;
      ref[num_objects] = .4; // low reflection
      hloid(-.1, 0, 1, .1);
      num_objects++;


  		double a,b,c ;
  		double H = Tan_half_angle;
  		a = Half_window_size / H ;
  		b = 0.5*window_width ;
  		c = 0.5*window_height ;
      double rgb[3];
  		double eyeTip[3];

      M3d_mat_mult_pt(eye, vm, eye); // put eye in eye space
  		for(double k = 0; k < 800; k+=1){
  			for(double l = 0; l < 800; l+=1){
          rgb[0] = world[0]; rgb[1] = world[1]; rgb[2] = world[2];

          eyeTip[0] = (k - b) / a;    // place the eyeTip relative to eye in eye space
  				eyeTip[1] = (l - c) / a;
          eyeTip[2] = 1;

  				ray (eye, eyeTip, rgb, 0) ;
  				G_rgb(rgb[0], rgb[1], rgb[2]) ;
  				G_point(k, l);
        }
      }

      M3d_mat_mult_pt(eye, vi, eye);// put eye back in world space

      G_rgb(0, 1, 0);
      G_circle(400 - 2.5, 400 - 2.5, 5);
      q = G_wait_key(); // wait for some direction on translation/rotation
      if(q == 65362){// up
        // printf("up\n");
        eye[0] += del*v[0]; // make eye go in direction of the vector
        eye[1] += del*v[1];
        eye[2] += del*v[2];
      }else if(q == 65363){
        angle += del*5;
        v[0] = sin(angle); // rotate the vector from original direction to new angle
        v[1] = 0;
        v[2] = cos(angle);
        // printf("right\n");
      }else if(q == 65364){
        eye[0] -= del*v[0];
        eye[1] -= del*v[1];
        eye[2] -= del*v[2];
        // printf("down\n");
      }else if(q == 65361){
        // printf("left\n");
        angle -= del*5;
        v[0] = sin(angle);
        v[1] = 0;
        v[2] = cos(angle);
      }
  }


}
