// --------------------------------------------------------------------------
// gMini,
// a minimal Glut/OpenGL app to extend                              
//
// Copyright(C) 2007-2009                
// Tamy Boubekeur
//                                                                            
// All rights reserved.                                                       
//                                                                            
// This program is free software; you can redistribute it and/or modify       
// it under the terms of the GNU General Public License as published by       
// the Free Software Foundation; either version 2 of the License, or          
// (at your option) any later version.                                        
//                                                                            
// This program is distributed in the hope that it will be useful,            
// but WITHOUT ANY WARRANTY; without even the implied warranty of             
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              
// GNU General Public License (http://www.gnu.org/licenses/gpl.txt)           
// for more details.                                                          
//                                                                          
// --------------------------------------------------------------------------

uniform float ambientRef;
uniform float diffuseRef;
uniform float specularRef;
uniform float shininess;

varying vec4 p;
varying vec3 n;

void main (void) {
    vec3 P = vec3 (gl_ModelViewMatrix * p); //Position du point à éclairer
    vec3 N = normalize (gl_NormalMatrix * n); //Normal en ce point
    vec3 V = normalize (-P); //Vecteur de vue
    
    vec4 Isa = gl_LightModel.ambient;
    vec4 Ka = gl_FrontMaterial.ambient;
    vec4 Ia = Isa * Ka;

    vec4 I = ambientRef * Ia;
  
    vec4 Isd = gl_LightSource[0].diffuse;
    vec4 Kd = gl_FrontMaterial.diffuse;
    vec3 vec_ps = vec3 (gl_LightSource[0].position.xyz - P);
    vec_ps = normalize(vec_ps);
    vec4 Id = Isd * Kd * dot(vec_ps, P);

    I += diffuseRef * Id;

    vec4 Iss = gl_LightSource[0].specular;
    vec4 Ks = gl_FrontMaterial.specular;

    vec3 R = vec3 ((2. * dot(vec_ps, P) * vec_ps) - P);
    R = normalize(R);

    vec4 Is = Iss * Ks * pow(dot(R, V), shininess);

    I += specularRef * Is;

    ////////////////////////////////////////////////
    //Eclairage de Phong à calculer en utilisant
    ////////////////////////////////////////////////
    // gl_LightSource[i].position.xyz Position de la lumière i
    // gl_LightSource[i].diffuse Couleur diffuse de la lumière i
    // gl_LightSource[i].specular Couleur speculaire de la lumière i
    // gl_FrontMaterial.diffuse Matériaux diffus de l'objet
    // gl_FrontMaterial.specular Matériaux speculaire de l'objet


    gl_FragColor = vec4 (I.xyz, 1);
}



