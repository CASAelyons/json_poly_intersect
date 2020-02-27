#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>
#include <time.h>
#include <argp.h>
#include <jansson.h>

static struct argp_option options[] = {
  {"output_geojson",   'j', "output_geojson_file",      0,   "Generate a geoJSON file" },
  {"geojson",          'g', "input_geojson_file",       0,   "Specify name of input geoJSON file"},
  { 0 }
};

#define MAX_NAME 1024

struct latLon {
  double lat;
  double lon;
};

struct arguments {
  char main_filename[MAX_NAME];
  char output_json_filename[MAX_NAME];
  char input_json_filename[MAX_NAME];
};

static error_t parse_opt(int key, char *arg, struct argp_state *state){
  struct arguments *arguments = state->input;

  switch(key){
    case 'j':
    if (arg[0] != '-')
      strncpy(arguments->output_json_filename,arg,MAX_NAME);
    else {
      printf("arguments should not begin with a - sign.  Exiting...\n");
      exit(0);
     }
    break;
  case 'g':
    if (arg[0] != '-')
      strncpy(arguments->input_json_filename,arg,MAX_NAME);
    else {
      printf("arguments should not begin with a - sign.  Exiting...\n");
      exit(0);
    }
    break;
  case ARGP_KEY_ARG:
    if(state->arg_num>=1)
      argp_usage(state);
    strncpy(arguments->main_filename,arg,MAX_NAME);
    break;
  case ARGP_KEY_END:
    if(state->arg_num<1)
      argp_usage(state);
    break;
  default:
    return ARGP_ERR_UNKNOWN;
  }
  return 0;
}

static struct argp argp = {options,parse_opt,"<input_json_file>","$Id: qpeExtract.c,v 1.0 2018-05-14 11:09:23 elyons Exp $"};

double min (double v1, double v2) {
  if (v1 < v2)
    return v1;
  else
    return v2;
}

double max (double v1, double v2) {
  if (v1 > v2)
    return v1;
  else
    return v2;
}

int pnpoly(int nvert, double *vertx, double *verty, double testx, double testy)
{
  /*
    Copyright (c) 1970-2003, Wm. Randolph Franklin
    Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

    Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimers.
    Redistributions in binary form must reproduce the above copyright notice in the documentation and/or other materials provided with the distribution.
    The name of W. Randolph Franklin may not be used to endorse or promote products derived from this Software without specific prior written permission. 

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
  */
  int i, j, c = 0;
  for (i = 0, j = nvert-1; i < nvert; j = i++) {
    if ( ((verty[i]>testy) != (verty[j]>testy)) &&
	 (testx < (vertx[j]-vertx[i]) * (testy-verty[i]) / (verty[j]-verty[i]) + vertx[i]) )
      c = !c;
  }
  return c;
}

int equal(double val1, double val2)
{
  if (fabs(val1 - val2) < .000000000001)
    return 1;
  else
    return 0;
}

int nearly_equal(double val1, double val2)
{
  if (fabs(val1 - val2) < .000000001)
    return 1;
  else
    return 0;
}


int getIntersection(struct latLon p1, struct latLon p2, struct latLon q1, struct latLon q2, struct latLon *intersection)
{
  double dlat1 = p2.lat - p1.lat;
  double dlon1 = p1.lon - p2.lon;
  double f1 = (dlat1 * p1.lon) + (dlon1 * p1.lat);
  //printf("dlat1: %f dlon1: %f f1: %f\n", dlat1, dlon1, f1);
  
  double dlat2 = q2.lat - q1.lat;
  double dlon2 = q1.lon - q2.lon;
  double f2 = (dlat2 * q1.lon) + (dlon2 * q1.lat);
  //printf("dlat2: %f dlon2: %f f2: %f\n", dlat2, dlon2, f2);
  
  double cp = (dlat1 * dlon2) - (dlat2 * dlon1);
  //printf("cp: %f ", cp);
  if (equal(cp, 0)) {
    //lines are parallel
    //printf("lines are parallel\n");
    return 0;
  }
  else {
    double ilon = ((dlon2 * f1) - (dlon1 * f2)) / cp;
    double ilat = ((dlat1 * f2) - (dlat2 * f1)) / cp;
    int valid1 = 0;
    int valid2 = 0;
    if (((min(p1.lon, p2.lon) < ilon) || equal(min(p1.lon, p2.lon), ilon)) &&
	((max(p1.lon, p2.lon) > ilon) || equal(max(p1.lon, p2.lon), ilon)) &&
	((min(p1.lat, p2.lat) < ilat) || equal(min(p1.lat, p2.lat), ilat)) &&
	((max(p1.lat, p2.lat) > ilat) || equal(min(p1.lat, p2.lat), ilat))) {
      valid1 = 1;
      //printf("valid1 ");
    }
    else
      valid1 = 0;

    if (((min(q1.lon, q2.lon) < ilon) || equal(min(q1.lon, q2.lon), ilon)) &&
	((max(q1.lon, q2.lon) > ilon) || equal(max(q1.lon, q2.lon), ilon)) &&
	((min(q1.lat, q2.lat) < ilat) || equal(min(q1.lat, q2.lat), ilat)) &&
	((max(q1.lat, q2.lat) > ilat) || equal(min(q1.lat, q2.lat), ilat))) {
      valid2 = 1;
      //printf("valid2 ");
    }
    else
      valid2 = 0;

    if ((valid1) && (valid2)) {
      intersection->lon = ilon;
      intersection->lat = ilat;
      //printf("intersect!\n");
      return 1;
    }
    else 
      return 0;
  }
}
      
int main(int argc, char* argv[])
{

  struct arguments arguments;
  arguments.output_json_filename[0] = '\0';
  arguments.input_json_filename[0] = '\0';
  arguments.main_filename[0] = '\0';
  
  argp_parse(&argp,argc,argv,0,0,&arguments);
  
  if ((arguments.input_json_filename[0] == '\0') || (arguments.main_filename[0] == '\0') || (arguments.output_json_filename[0] == '\0')) {
      printf("Please specify arguments correctly... exiting\n");
      exit(0);
  }
  
  //main file read
  json_t *json;
  json_error_t error;
  json_t *features;
  size_t num_features;
  int k;
  
  //Ok, with that, we're off...
  json = json_load_file(arguments.main_filename, 0, &error);
  if(!json) {
    printf("JSON ERROR: %s\n", error.text);
    exit(-1);
  }
  
  if (!json_is_object(json)){
    printf("json type: %d\n", json_typeof(json));
    printf("geoJSON file does not contain an object.  Please try another geoJSON file. \n");
    exit(0);
  }
  
  features = json_object_get(json, "features");
  if (features == NULL) {
    printf("geoJSON file does not contain the necessary \"features\" field.  Exiting...\n");
    exit(0);
  }
  if (!json_is_array(features)) {
    printf("\"features\" field is not an array.  Exiting...");
    exit(0);
  }

  num_features = json_array_size(features);
  if (num_features < 1) {
    printf("no features present in this geoJSON.  Exiting...");
    exit(0);
  }
  
  //Appears we have a good main geoJSON file.  Lets open the input geoJSON file...
  json_t *injson;
  json_t *infeatures;
  size_t num_infeatures;
  int l;
  injson = json_load_file(arguments.input_json_filename, 0, &error);
  if(!injson) {
    printf("INJSON ERROR: %s\n", error.text);
    exit(-1);
  }

  if (!json_is_object(injson)){
    printf("json type: %d\n", json_typeof(injson));
    printf("input_geoJSON file does not contain an object.  Please try another geoJSON file. \n");
    exit(0);
  }

  infeatures = json_object_get(injson, "features");
  if (infeatures == NULL) {
    printf("input_geoJSON file does not contain the necessary \"features\" field.  Exiting...\n");
    exit(0);
  }
  if (!json_is_array(infeatures)) {
    printf("\"features\" field is not an array.  Exiting...");
    exit(0);
  }

  num_infeatures = json_array_size(infeatures);
  if (num_infeatures < 1) {
    printf("no features present in this input_geoJSON.  Exiting...");
    exit(0);
  }

  //ok that one is good too....  Lets go
  
  for (k=0; k<num_features; k++) {
    json_t *this_feature = json_array_get(features, k);
    
    if (json_is_object(this_feature)) {
      json_t *this_feature_type = json_object_get(this_feature, "type");
      if (json_is_string(this_feature_type)) { 
	if (strcmp(json_string_value(this_feature_type), "Feature") == 0) {
	  json_t *this_feature_geometry = json_object_get(this_feature, "geometry");
	  if (json_is_object(this_feature_geometry)) {
	    json_t *this_feature_geometry_type = json_object_get(this_feature_geometry, "type");
	    if (json_is_string(this_feature_geometry_type)) {
	      int is_poly = 0;
	      int is_multipoly = 0;
	      if (strcmp(json_string_value(this_feature_geometry_type), "Polygon") == 0)
		is_poly = 1;
	      else if (strcmp(json_string_value(this_feature_geometry_type), "MultiPolygon") == 0)
		is_multipoly = 1;
	      
	      if ((is_poly) || (is_multipoly)) {
		json_t *this_feature_geometry_coordinates = json_object_get(this_feature_geometry, "coordinates");
		if (json_is_array(this_feature_geometry_coordinates)) {
		  size_t num_polygons = json_array_size(this_feature_geometry_coordinates);
		  int p,q;
		  
		  if (is_poly) {	    
		    for (p=0; p<num_polygons; p++) {
		      json_t *this_polygon = json_array_get(this_feature_geometry_coordinates, p);
		      if (json_is_array(this_polygon)) {
			size_t num_coords = json_array_size(this_polygon);
			int c;
			double *latpoints, *lonpoints;
			size_t goodvals = 0;
			latpoints = malloc(num_coords * sizeof(double));
			lonpoints = malloc(num_coords * sizeof(double));
			for (c=0; c<num_coords; c++) {
			  json_t *llpair = json_array_get(this_polygon, c);
			  if (json_is_array(llpair)) {
			    if (json_array_size(llpair) == 2) {
			      if ((json_is_real(json_array_get(llpair, 0))) && (json_is_real(json_array_get(llpair, 1)))) {
				lonpoints[c] = json_real_value(json_array_get(llpair, 0));
				latpoints[c] = json_real_value(json_array_get(llpair, 1));
				goodvals++;
			      }
			    }
			  }
			}
			if (goodvals == num_coords) {
			  //main file is good and this features points are loaded into latpoints/lonpoints array.
			  
			  for (l=0; l<num_infeatures; l++) {
			    //printf("k: %d l: %d\n", k, l);
			    json_t *this_infeature = json_array_get(infeatures, l);
			    
			    if (json_is_object(this_infeature)) {
			      json_t *this_infeature_type = json_object_get(this_infeature, "type");
			      if (json_is_string(this_infeature_type)) {
				if (strcmp(json_string_value(this_infeature_type), "Feature") == 0) {
				  json_t *this_infeature_properties = json_object_get(this_infeature, "properties");
				  //const char *key;
				  //json_t *value;
				  //json_object_foreach(this_infeature_properties, key, value) { printf("key %s val %s\n", key, json_string_value(value));}
				  json_t *this_infeature_name = json_object_get(this_infeature_properties, "NAME");
				  printf("This feature name: %s\n", json_string_value(this_infeature_name));				    
				  json_t *this_infeature_geometry = json_object_get(this_infeature, "geometry");
				  if (json_is_object(this_infeature_geometry)) {
				    json_t *this_infeature_geometry_type = json_object_get(this_infeature_geometry, "type");
				    if (json_is_string(this_infeature_geometry_type)) {
				      int in_is_poly = 0;
				      int in_is_multipoly = 0;
				      if (strcmp(json_string_value(this_infeature_geometry_type), "Polygon") == 0)
					in_is_poly = 1;
				      else if (strcmp(json_string_value(this_infeature_geometry_type), "MultiPolygon") == 0)
					in_is_multipoly = 1;

				      if ((in_is_poly) || (in_is_multipoly)) {
					json_t *this_infeature_geometry_coordinates = json_object_get(this_infeature_geometry, "coordinates");
					if (json_is_array(this_infeature_geometry_coordinates)) {
					  size_t num_inpolygons = json_array_size(this_infeature_geometry_coordinates);
					  int in_p,in_q;
					  if (in_is_poly) {
					    for (in_p=0; in_p<num_inpolygons; in_p++) {
					      json_t *ints_arr_json = json_array();
					      json_t *this_inpolygon = json_array_get(this_infeature_geometry_coordinates, in_p);
					      if (json_is_array(this_inpolygon)) {
						size_t num_incoords = json_array_size(this_inpolygon);
						int in_c;
						double *inlatpoints, *inlonpoints;
						size_t ingoodvals = 0;
						inlatpoints = malloc(num_incoords * sizeof(double));
						inlonpoints = malloc(num_incoords * sizeof(double));
						for (in_c=0; in_c<num_incoords; in_c++) {
						  json_t *llpair = json_array_get(this_inpolygon, in_c);
						  if (json_is_array(llpair)) {
						    if (json_array_size(llpair) == 2) {
						      if ((json_is_real(json_array_get(llpair, 0))) && (json_is_real(json_array_get(llpair, 1)))) {
							inlonpoints[in_c] = json_real_value(json_array_get(llpair, 0));
							inlatpoints[in_c] = json_real_value(json_array_get(llpair, 1));
							ingoodvals++;
						      }
						    }
						  }
						}

						if (ingoodvals == num_incoords) {
						  int broach = -1;
						  int *cross_in_point_arr;
						  int *cross_out_point_arr;
						  int *main_cross_in_point_arr;
						  int *main_cross_out_point_arr;
						  int tmp;
						  
						  cross_in_point_arr = (int*)malloc(sizeof(int));
						  cross_out_point_arr = (int*)malloc(sizeof(int));
						  cross_in_point_arr[0] = -1;
						  cross_out_point_arr[0] = -1;
						  main_cross_in_point_arr = (int*)malloc(sizeof(int));
						  main_cross_out_point_arr =(int*)malloc(sizeof(int));
						  main_cross_in_point_arr[0] = -1;
						  main_cross_out_point_arr[0] = -1;
						  int num_cross_ins = 0;
						  int num_cross_outs = 0;
						  
						  int process_data = 1;
						  int end_c = num_incoords + 1;
						  in_c = 0;

						  //simple check of the first point...
						  //tells us if the cross in points will be known in advance of the cross out points
						  int initial;
						  if (pnpoly(num_coords, lonpoints, latpoints, inlonpoints[0], inlatpoints[0]))
						    initial = 1;
						  else
						    initial = 0;
						  
						  while (in_c != end_c) {
						    struct latLon interPoint;
						    struct latLon inPoint1;
						    struct latLon inPoint2;
						    struct latLon mainPoint1;
						    struct latLon mainPoint2;
						    struct latLon *intersection_array;
						    intersection_array = (struct latLon*)malloc(sizeof(struct latLon));
						    int edgecount = 0;
						    while (process_data) {
						      int inside = pnpoly(num_coords, lonpoints, latpoints, inlonpoints[in_c], inlatpoints[in_c]);
						      if (inside == 1) {
							if (broach == 0) {
							  //start outside, cross in
							  //first add the intersection
							  inPoint1.lat = inlatpoints[in_c-1];
							  inPoint1.lon = inlonpoints[in_c-1];
							  inPoint2.lat = inlatpoints[in_c];
							  inPoint2.lon = inlonpoints[in_c];
							  for (c=0; c<num_coords - 1; c++) {
							    mainPoint1.lat = latpoints[c];
							    mainPoint1.lon = lonpoints[c];
							    mainPoint2.lat = latpoints[c+1];
							    mainPoint2.lon = lonpoints[c+1];
							    if (getIntersection(mainPoint1, mainPoint2, inPoint1, inPoint2, &interPoint)) {
							      if (pnpoly(num_incoords, inlonpoints, inlatpoints, lonpoints[c], latpoints[c])) {
								tmp = c;
							      }
							      else if (pnpoly(num_incoords, inlonpoints, inlatpoints, lonpoints[c+1], latpoints[c+1])) {
								tmp = c+1;
							      }
							      else {
								printf("points overlap... exiting\n");
								exit(-1);
							      }
							      cross_in_point_arr[num_cross_ins] = in_c;
							      main_cross_in_point_arr[num_cross_ins] = tmp;
							      
							      ++num_cross_ins;
							      cross_in_point_arr = (int*)realloc(cross_in_point_arr, (num_cross_ins+1)*sizeof(int));
							      main_cross_in_point_arr = (int*)realloc(main_cross_in_point_arr, (num_cross_ins+1)*sizeof(int));
							
							      intersection_array[edgecount].lat = interPoint.lat;
							      intersection_array[edgecount].lon = interPoint.lon;
							      ++edgecount;
							      intersection_array = (struct latLon*)realloc(intersection_array, (edgecount+1)*sizeof(struct latLon));
							      c = num_coords - 1;
							    }
							  }
							  //now add the current point which is inside
							  intersection_array[edgecount].lat = inlatpoints[in_c];
							  intersection_array[edgecount].lon = inlonpoints[in_c];
							  ++edgecount;
							  intersection_array = (struct latLon*)realloc(intersection_array, (edgecount+1)*sizeof(struct latLon));
							}
							else {
							  //inside... add the point
							  intersection_array[edgecount].lat = inlatpoints[in_c];
							  intersection_array[edgecount].lon = inlonpoints[in_c];
							  ++edgecount;
							  intersection_array = (struct latLon*)realloc(intersection_array, (edgecount+1)*sizeof(struct latLon));
							}
							//denote that this, the most recent point, is inside the polygon
							broach = 1;
						      }
						      else {
							//this point is not inside the main polygon
							if (broach == 1) {
							  //started inside, crossed out
							  int direction = 0;
							  int in_c_direction = 0;
							  int start_point;
							  //add the intersection
							  inPoint1.lat = inlatpoints[in_c-1];
							  inPoint1.lon = inlonpoints[in_c-1];
							  inPoint2.lat = inlatpoints[in_c];
							  inPoint2.lon = inlonpoints[in_c];
							  for (c=0; c<num_coords - 1; c++) {
							    mainPoint1.lat = latpoints[c];
							    mainPoint1.lon = lonpoints[c];
							    mainPoint2.lat = latpoints[c+1];
							    mainPoint2.lon = lonpoints[c+1];
							    if (getIntersection(mainPoint1, mainPoint2, inPoint1, inPoint2, &interPoint)) {
							      cross_out_point_arr[num_cross_outs] = in_c;
							      if (pnpoly(num_incoords, inlonpoints, inlatpoints, lonpoints[c], latpoints[c])) {
								tmp = c;
								if (tmp == 0)
								  start_point = num_coords - 1;
								else
								  start_point = tmp - 1;
								direction = -1;
							      }
							      else if (pnpoly(num_incoords, inlonpoints, inlatpoints, lonpoints[c+1], latpoints[c+1])) {
								tmp = c+1;
								if (tmp == num_coords - 1)
								  start_point = 0;
								else
								  start_point = tmp + 1;
								direction = 1;
							      }
							      else {
								printf("logic error on cross in points.  exiting...\n");
								exit(-1);
							      }
							      main_cross_out_point_arr[num_cross_outs] = tmp;
							      
							      ++num_cross_outs;
							      cross_out_point_arr = (int*)realloc(cross_out_point_arr, (num_cross_outs+1)*sizeof(int));
							      main_cross_out_point_arr = (int*)realloc(main_cross_out_point_arr, (num_cross_outs+1)*sizeof(int));

							      intersection_array[edgecount].lat = interPoint.lat;
							      intersection_array[edgecount].lon = interPoint.lon;
							      ++edgecount;
							      intersection_array = (struct latLon*)realloc(intersection_array, (edgecount+1)*sizeof(struct latLon));
							      c = num_coords - 1;
							    }
							  }
							  //now add the main point and head in the direction that continues inside the input poly
					   
							  tmp = main_cross_out_point_arr[num_cross_outs - 1];
							  intersection_array[edgecount].lat = latpoints[tmp];
							  intersection_array[edgecount].lon = lonpoints[tmp];
							  ++edgecount;
							  intersection_array = (struct latLon*)realloc(intersection_array, (edgecount+1)*sizeof(struct latLon));

							  int old_start_point;
							  int in_start_point = -1;
							  //two termination conditions:
							  //1, that start point gets back to the original main cross in point
							  //2, that in_start_point gets to zero
							  
							  int continue_loop = 1;
							  while(continue_loop) {
							  
							    //proceed down the main polygon, checking for intersections
							    while (pnpoly(num_incoords, inlonpoints, inlatpoints, lonpoints[start_point], latpoints[start_point])) {
							      intersection_array[edgecount].lat = latpoints[start_point];
							      intersection_array[edgecount].lon = lonpoints[start_point];
							      ++edgecount;
							      intersection_array = (struct latLon*)realloc(intersection_array, (edgecount+1)*sizeof(struct latLon));
							      old_start_point = start_point;
							      start_point = start_point + direction;
							      if (start_point == num_coords)
								start_point = 0;
							      if (start_point == -1)
								start_point = num_coords - 1;
							    }
							    //you've come to an intersection or the end
							    //case a, the end...
							    if (old_start_point == main_cross_in_point_arr[0]) {
							      
							      in_c = num_incoords - 1;
							      continue_loop = 0;
							    }
							    else {
							      //add the point that you crossed out to the cross out array
							      main_cross_in_point_arr[num_cross_ins] = old_start_point;
							     
							      //add the intersection to the output polygon
							      mainPoint1.lat = latpoints[start_point];
							      mainPoint1.lon = lonpoints[start_point];
							      mainPoint2.lat = latpoints[main_cross_in_point_arr[num_cross_ins]];
							      mainPoint2.lon = lonpoints[main_cross_in_point_arr[num_cross_ins]];
							      for (c=0; c<num_incoords - 1; c++) {
								inPoint1.lat = inlatpoints[c];
								inPoint1.lon = inlonpoints[c];
								inPoint2.lat = inlatpoints[c+1];
								inPoint2.lon = inlonpoints[c+1];
								if (getIntersection(inPoint1, inPoint2, mainPoint1, mainPoint2, &interPoint)) {
								  if (pnpoly(num_coords, lonpoints, latpoints, inlonpoints[c], inlatpoints[c])) {
								    cross_in_point_arr[num_cross_ins] = c;
								    in_c_direction = -1;
								  }
								  else if (pnpoly(num_coords, lonpoints, latpoints, inlonpoints[c+1], inlatpoints[c+1])) {
								    cross_in_point_arr[num_cross_ins] = c+1;
								    in_c_direction = 1;
								  }
								  else {
								    printf("logic error on cross in points.  exiting...\n");
								    exit(-1);
								  }
								  ++num_cross_ins;
								  cross_in_point_arr = (int*)realloc(cross_in_point_arr, (num_cross_ins+1)*sizeof(int));
								  main_cross_in_point_arr = (int*)realloc(main_cross_in_point_arr, (num_cross_ins+1)*sizeof(int));
								  
								  intersection_array[edgecount].lat = interPoint.lat;
								  intersection_array[edgecount].lon = interPoint.lon;
								  ++edgecount;
								  intersection_array = (struct latLon*)realloc(intersection_array, (edgecount+1)*sizeof(struct latLon));
								  c = num_incoords - 1;
								}
							      }
							      
							      //repeat the searching on the input poly
							      in_start_point = cross_in_point_arr[num_cross_ins - 1];
							      
							     
							      if ((initial) && (in_start_point == 0)) {}
							      else {
								while (pnpoly(num_coords, lonpoints, latpoints, inlonpoints[in_start_point], inlatpoints[in_start_point])) {
								  if ((initial) && (in_start_point == 0))
								    break;
								  intersection_array[edgecount].lat = inlatpoints[in_start_point];
								  intersection_array[edgecount].lon = inlonpoints[in_start_point];
								  ++edgecount;
								  intersection_array = (struct latLon*)realloc(intersection_array, (edgecount+1)*sizeof(struct latLon));
								  old_start_point = in_start_point;
								  in_start_point = in_start_point + in_c_direction;
								  if (in_start_point == num_incoords)
								    in_start_point = 0;
								  if (in_start_point == -1)
								    in_start_point = num_incoords - 1;
								  
								}
							      }
							      if ((initial) && (in_start_point == 0)) {
								in_c = num_incoords - 1;
							      }
							      else {
								//back to an intersection
							      
								inPoint1.lat = inlatpoints[old_start_point];
								inPoint1.lon = inlonpoints[old_start_point];
								inPoint2.lat = inlatpoints[in_start_point];
								inPoint2.lon = inlonpoints[in_start_point];
								for (c=0; c<num_coords - 1; c++) {
								  mainPoint1.lat = latpoints[c];
								  mainPoint1.lon = lonpoints[c];
								  mainPoint2.lat = latpoints[c+1];
								  mainPoint2.lon = lonpoints[c+1];
								  if (getIntersection(mainPoint1, mainPoint2, inPoint1, inPoint2, &interPoint)) {
								    cross_out_point_arr[num_cross_outs] = in_start_point;
								    if (pnpoly(num_incoords, inlonpoints, inlatpoints, lonpoints[c], latpoints[c])) {
								      tmp = c;
								      if (tmp == 0)
									start_point = num_coords - 1;
								      else
									start_point = tmp - 1;
								      direction = -1;
								    }
								    else if (pnpoly(num_incoords, inlonpoints, inlatpoints, lonpoints[c+1], latpoints[c+1])) {
								      tmp = c+1;
								      if (tmp == num_coords - 1)
									start_point = 0;
								      else
									start_point = tmp + 1;
								      direction = 1;
								    }
								    else {
								      printf("logic error on cross in points.  exiting...\n");
								      exit(-1);
								    }
								    main_cross_out_point_arr[num_cross_outs] = tmp;
								    
								    ++num_cross_outs;
								    cross_out_point_arr = (int*)realloc(cross_out_point_arr, (num_cross_outs+1)*sizeof(int));
								    main_cross_out_point_arr = (int*)realloc(main_cross_out_point_arr, (num_cross_outs+1)*sizeof(int));
								    
								    intersection_array[edgecount].lat = interPoint.lat;
								    intersection_array[edgecount].lon = interPoint.lon;
								    ++edgecount;
								    intersection_array = (struct latLon*)realloc(intersection_array, (edgecount+1)*sizeof(struct latLon));
								    c = num_coords - 1;
								  }
								}
							      }
							    }
							    if ((old_start_point == main_cross_in_point_arr[0]) || ((initial) && in_start_point == 0)) {
							      
							      continue_loop = 0;
							    }
							    else {
							      //loop the while array
							      
							    }
							    
							  }

							  process_data = 0;
							  //ok we've got a polygon here ready to go
							}
							else {
							  //last point was outside, still outside.. do nothing
							}
							//denote that this point was outside
							broach = 0; 
						      }
						      //increment the counter
						      ++in_c;

						      //check if we're done or need to wrap back around
						      if ((in_c == num_incoords) && (end_c == num_incoords + 1)) {
							process_data = 0;
							++in_c;
						      }
						      else if (in_c == end_c)
							process_data = 0;
						      else if (in_c == num_incoords) 
							in_c = 0;
						    }
						  
						    //printf("at the end of this poly we now have %d points\n", edgecount);
						  
						  //at this point you have all the points for the intersecting polygon in the intersection array
						  //they are unordered and may contain duplicates however.... let's remove any duplicates
						    
						    if (edgecount > 2) { 
						      struct latLon *ints_arr;
						      int final_points = 0;
						      ints_arr = (struct latLon*)malloc(sizeof(struct latLon));
						      for (c = 0; c < edgecount-1; c++) {
							int duplicate_found = 0;
							int tmp_c;
							for (tmp_c = c+1; tmp_c < edgecount; tmp_c++) {
							  if ((nearly_equal(intersection_array[c].lat, intersection_array[tmp_c].lat)) && (nearly_equal(intersection_array[c].lon, intersection_array[tmp_c].lon))){
							    duplicate_found = 1;
							    tmp_c = edgecount;
							  }
							}
							if (duplicate_found == 0) {
							  ints_arr[final_points] = intersection_array[c];
							  final_points++;
							  ints_arr = (struct latLon*)realloc(ints_arr, (final_points+1)*sizeof(struct latLon));
							}
						      }
						      
						      //for (c = 0; c < final_points; c++)
						      //  printf("c: %d ints_arr_lat: %f ints_arr_lon: %f\n", c, ints_arr[c].lat, ints_arr[c].lon);
						    
						      json_t *ints_poly_coord_array = json_array();
						      for (c = 0; c < final_points; c++){
							json_t *ints_coord_array = json_array();
							json_t *coord_lon = json_real(ints_arr[c].lon);
							json_t *coord_lat = json_real(ints_arr[c].lat);
							json_array_append(ints_coord_array, coord_lon);
							json_array_append(ints_coord_array, coord_lat);
							//printf("c: %d ints_arr_lat: %d ints_arr_lon: %d\n", c, ints_arr[c].lat, ints_arr[c].lon);
							//printf("ints_coord_array 0: %f ints_coord_array 1: %f \n", json_real_value(json_array_get(ints_coord_array, 0)), json_real_value(json_array_get(ints_coord_array, 1)));
							//json_array_set(ints_arr_json, c, ints_coord_array); 
							json_array_append(ints_poly_coord_array, ints_coord_array);
						      }
						      
						      //append the last point, same as the first point
						      json_t *ints_coord_array = json_array();
						      json_t *coord_lon = json_real(ints_arr[0].lon);
						      json_t *coord_lat = json_real(ints_arr[0].lat);
						      json_array_append(ints_coord_array, coord_lon);
						      json_array_append(ints_coord_array, coord_lat);
						      json_array_append(ints_poly_coord_array, ints_coord_array);
						      
						      //json_array_set(ints_arr_json, i, ints_poly_coord_array);
						      
						      //json_array_clear(this_infeature_geometry_coordinates);
						      //json_array_set(this_inpolygon, 
						      //json_array_set(this_infeature_geometry_coordinates, in_p, ints_arr_json);
						      //json_array_set(this_infeature_geometry_coordinates, in_p, ints_arr_json);
						      //json_object_set(this_infeature_geometry, "coordinates", ints_arr_json);
						      
						      json_array_append(ints_arr_json, ints_poly_coord_array);
						      free(intersection_array);
						      process_data = 1;
						      //json_array_set(this_infeature_geometry_coordinates, in_p, ints_arr_json);
						      
						      //json_array_set(this_infeature_geometry_coordinates, in_p, ints_poly_coord_array);
						    }
						  }
						  size_t num_polys = json_array_size(ints_arr_json);
						  printf("num_polys: %d\n", (int)num_polys);
						  json_array_set(this_infeature_geometry_coordinates, in_p, ints_arr_json);
						  //if (num_polys > 1) {
						  json_t *jsonmp = json_string("MultiPolygon");
						  json_object_set(this_infeature_geometry, "type", jsonmp);
						    //}
						  free(cross_in_point_arr);
						  free(cross_out_point_arr);
						  free(main_cross_in_point_arr);
						  free(main_cross_out_point_arr);
						}
						//
					      }
					    }
					  }
					}
					//json_object_del(this_infeature_geometry, "coordinates");
					//json_object_set(this_infeature_geometry, "coordinates", this_infeature_geometry_coordinates);
				      }
				    }
				  }
				}
			      }
			    }
			  }
			}
		      }
		    }
		  }
		  else if (is_multipoly) {
		    /*
		      for (p=0; p<num_polygons; p++) {
		      json_t *this_multipolygon = json_array_get(this_feature_geometry_coordinates, p);
		      if (json_is_array(this_multipolygon)) {
		      size_t num_polys = json_array_size(this_multipolygon);
		      for (q=0; q<num_polys; q++) {
		      json_t *this_multipoly = json_array_get(this_multipolygon, q);
		      if (json_is_array(this_multipoly)) {
		      size_t num_coords = json_array_size(this_multipoly);
		      int c;
		      double *latpoints, *lonpoints;
		      size_t goodvals = 0;
		      latpoints = malloc(num_coords * sizeof(double));
		      lonpoints = malloc(num_coords * sizeof(double));
		      for (c=0; c<num_coords; c++) {
		      json_t *llpair = json_array_get(this_multipoly, c);
		      if (json_is_array(llpair)) {
		      if (json_array_size(llpair) == 2) {
		      if ((json_is_real(json_array_get(llpair, 0))) && (json_is_real(json_array_get(llpair, 1)))) {
		      lonpoints[c] = json_real_value(json_array_get(llpair, 0));
		      latpoints[c] = json_real_value(json_array_get(llpair, 1));
		      goodvals++;
		      }
		      }
		      }
		      }
		      if (goodvals == num_coords) {
		      for (i=0; i<num_lats; i++) {
		      for (j=0; j<num_lons; j++) {
		      int inside = pnpoly(num_coords, lonpoints, latpoints, x_in[j], y_in[i]);
		      if (inside == 1) {
		      inpoly[i][j] = 1;
		      qpe_out[(num_lats - i) - 1][j] = qpe_in[(num_lats - i) - 1][j];
		      if ((arguments.averaging) && (qpe_in[(num_lats - i) - 1][j] > 1)) {		      num_inside_points++;
		      total = total + qpe_in[(num_lats - i) - 1][j];
		      }
		      }
		      }
		      }
		      }
		      }
		      }
		      }
		      }
		    */
		  }
		  
		  /*
		    if(arguments.property) {
		      int setaccumstatus = json_object_set(this_feature_properties, varname, json_real((double)polygon_average_IN)); 
		      int settimestatus = json_object_set(this_feature_properties, "Timestamp", json_string(starttime));
		      if ((setaccumstatus == -1) || (settimestatus == -1)) {
			printf("GeoJSON property append failed!\n");
		      }
		    */
		}
	      }
	    }
	  } 
	}
      }
    }
  }

  int dumpstatus = json_dump_file(injson, arguments.output_json_filename, 0);
  if (dumpstatus == -1) {
    printf("GeoJSON file dump failed!\n");
  }
 
return 0;
}
