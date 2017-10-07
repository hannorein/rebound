shader_code = """
<script id="orbit_shader-vs" type="x-shader/x-vertex">
    uniform vec3 focus;
    uniform vec3 aef;
    uniform vec3 omegaOmegainc;
    attribute float lintwopi;
    varying float lin;
    uniform mat4 mvp;
    const float M_PI = 3.14159265359;
    void main() {
       float a = aef.x;
       float e = aef.y;
       float f = aef.z+lintwopi;
       lin = lintwopi/(M_PI*2.);
       if (e>1.){
           float theta_max = acos(-1./e);
           f = 0.0001-theta_max+1.9998*lin*theta_max;
           lin = sqrt(min(0.5,lin));
       }
       float omega = omegaOmegainc.x;
       float Omega = omegaOmegainc.y;
       float inc = omegaOmegainc.z;
       float r = a*(1.-e*e)/(1. + e*cos(f));
       float cO = cos(Omega);
       float sO = sin(Omega);
       float co = cos(omega);
       float so = sin(omega);
       float cf = cos(f);
       float sf = sin(f);
       float ci = cos(inc);
       float si = sin(inc);
       vec3 pos = vec3(r*(cO*(co*cf-so*sf) - sO*(so*cf+co*sf)*ci),r*(sO*(co*cf-so*sf) + cO*(so*cf+co*sf)*ci),+ r*(so*cf+co*sf)*si);
       gl_Position = mvp*(vec4(focus+pos, 1.0));
    }
</script>
<script id="orbit_shader-fs" type="x-shader/x-fragment">
    precision mediump float;
    varying float lin;
    void main() {
      float fog = max(max(0.,-1.+2.*gl_FragCoord.z),max(0.,1.-2.*gl_FragCoord.z));
      gl_FragColor = vec4(1.,1.,1.,sqrt(lin)*(1.-fog));
    }
</script>
<script id="point_shader-vs" type="x-shader/x-vertex">
    attribute vec3 vp;
    uniform mat4 mvp;
    //uniform vec4 vc;
    //varying vec4 color;
    void main() {
      gl_PointSize = 15.0;
      gl_Position = mvp*vec4(vp, 1.0);
      //color = vc;
    }
</script>
<script id="point_shader-fs" type="x-shader/x-fragment">
    precision mediump float;
    //varying vec4 color;
    void main() {
      vec2 rel = gl_PointCoord.st;
      rel.s -=0.5;
      rel.t -=0.5;
      if (length(rel)>0.25){
         gl_FragColor = vec4(0.,0.,0.,0.); 
      }else{
         vec4 cmod = vec4(1.,1.,1.,1.);
         float fog = max(max(0.,-1.+2.*gl_FragCoord.z),max(0.,1.-2.*gl_FragCoord.z));
         cmod.a*= (1.-fog)*min(1.,1.-4.*(length(rel)/0.25-0.75));
         gl_FragColor = cmod;
      }
    }
</script>
"""
js_code = """
<script>
function compileShader(glr, shaderSource, shaderType) {
  // Create the shader object
  var shader = glr.createShader(shaderType);
 
  // Set the shader source code.
  glr.shaderSource(shader, shaderSource);
 
  // Compile the shader
  glr.compileShader(shader);
 
  // Check if it compiled
  var success = glr.getShaderParameter(shader, glr.COMPILE_STATUS);
  if (!success) {
    // Something went wrong during compilation; get the error
    throw "could not compile shader:" + glr.getShaderInfoLog(shader);
  }
 
  return shader;
}
function createShaderFromScript(glr, scriptId, opt_shaderType) {
  // look up the script tag by id.
  var shaderScript = document.getElementById(scriptId);
  if (!shaderScript) {
    throw("*** Error: unknown script element" + scriptId);
  }
 
  // extract the contents of the script tag.
  var shaderSource = shaderScript.text;
 
  // If we didn't pass in a type, use the 'type' from
  // the script tag.
  if (!opt_shaderType) {
    if (shaderScript.type == "x-shader/x-vertex") {
      opt_shaderType = glr.VERTEX_SHADER;
    } else if (shaderScript.type == "x-shader/x-fragment") {
      opt_shaderType = glr.FRAGMENT_SHADER;
    } else if (!opt_shaderType) {
      throw("*** Error: shader type not set");
    }
  }
 
  return compileShader(glr, shaderSource, opt_shaderType);
};
function createProgramFromScripts( glr, vertexShaderId, fragmentShaderId) {
  var vertexShader = createShaderFromScript(glr, vertexShaderId, glr.VERTEX_SHADER);
  var fragmentShader = createShaderFromScript(glr, fragmentShaderId, glr.FRAGMENT_SHADER);
  var program = glr.createProgram();
 
  // attach the shaders.
  glr.attachShader(program, vertexShader);
  glr.attachShader(program, fragmentShader);
 
  // link the program.
  glr.linkProgram(program);
 
  // Check if it linked.
  var success = glr.getProgramParameter(program, glr.LINK_STATUS);
  if (!success) {
      // something went wrong with the link
      throw ("program filed to link:" + glr.getProgramInfoLog (program));
  }
 
  return program;
}
function quat2mat(A,mat){
    var xx = A.x*A.x; var xy = A.x*A.y; var xz = A.x*A.z;
    var xw = A.x*A.w; var yy = A.y*A.y; var yz = A.y*A.z;
    var yw = A.y*A.w; var zz = A.z*A.z; var zw = A.z*A.w;
    mat[0] = 1.-2.*(yy+zz);
    mat[1] =    2.*(xy-zw);
    mat[2] =    2.*(xz+yw);
    mat[4] =    2.*(xy+zw);
    mat[5] = 1.-2.*(xx+zz);
    mat[6] =    2.*(yz-xw);
    mat[8] =    2.*(xz-yw);
    mat[9] =    2.*(yz+xw);
    mat[10]= 1.-2.*(xx+yy);
    mat[3] = mat[7] = mat[11] = mat[12] = mat[13] = mat[14] = 0.; mat[15]= 1.;
}
function multvec(A, B, vecr){
    var mat = [0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.];
    quat2mat(A,mat);
    vecr[0] = mat[0]*B[0] + mat[1]*B[1] + mat[2]*B[2];
    vecr[1] = mat[4]*B[0] + mat[5]*B[1] + mat[6]*B[2];
    vecr[2] = mat[8]*B[0] + mat[9]*B[1] + mat[10]*B[2];
}
function mattransp(mat){
    var matt = [
        mat[0], mat[4], mat[8], mat[12],
        mat[1], mat[5], mat[9], mat[13],
        mat[2], mat[6], mat[10], mat[14],
        mat[3], mat[7], mat[11], mat[15]];
    return matt;
}
function conjugate(quat){
    var cquat = {x:-quat.x, y:-quat.y, z:-quat.z, w:quat.w};
    return cquat;
}
function mult(A, B){
    var mquat = {   x: A.w*B.x + A.x*B.w + A.y*B.z - A.z*B.y,
                    y: A.w*B.y - A.x*B.z + A.y*B.w + A.z*B.x,
                    z: A.w*B.z + A.x*B.y - A.y*B.x + A.z*B.w,
                    w: A.w*B.w - A.x*B.x - A.y*B.y - A.z*B.z};
    return mquat;
}

function normalize(quat){
    var L = Math.sqrt(quat.x*quat.x + quat.y*quat.y + quat.z*quat.z + quat.w*quat.w);
    var nquat = {x:quat.x/L, y:quat.y/L, z:quat.z/L, w:quat.w/L};
    return nquat;
}
function matortho(mat, l, r, b, t, n, f){
    mat[0] = 2./(r-l); mat[1] = 0.; mat[2] = 0.; mat[3] = -(r+l)/(r-l);
    mat[4] = 0.; mat[5] = 2./(t-b); mat[6] = 0.; mat[7] = -(t+b)/(t-b);
    mat[8] = 0.; mat[9] = 0.; mat[10] = -2./(f-n); mat[11] = -(f+n)/(f-n);
    mat[12] = 0.; mat[13] = 0.; mat[14] = 0.; mat[15] = 1.;
}
function matmult(A,B,C){
    for(i=0;i<4;i++){
    for(j=0;j<4;j++){
        C[i+4*j] = 0.;
    for(k=0;k<4;k++){
        C[i+4*j] += A[k+4*j]*B[i+4*k];
    }}}
}
function startGL(reboundView) {
    var canvas = document.getElementById("reboundcanvas-"+reboundView.cid);
    if (!canvas){
        reboundView.startCount = reboundView.startCount+1;
        if (reboundView.startCount>1000){
            console.log("Cannot find element.");
        }else{
            setTimeout(function(){ startGL(reboundView); }, 10);
        }
        return;
    }
    var rect = canvas.getBoundingClientRect()
    reboundView.ratio = rect.width/rect.height;
    reboundView.view = normalize({x:reboundView.orientation[0], y:reboundView.orientation[1], z:reboundView.orientation[2], w:reboundView.orientation[3]});

    canvas.addEventListener('mousedown', function() {
        reboundView.mouseDown=1;
        }, false);
    canvas.addEventListener('mouseup', function() {
        reboundView.mouseDown=0;
        }, false);
    canvas.addEventListener('mouseleave', function() {
        reboundView.mouseDown=0;
        }, false);

    canvas.addEventListener('mousemove', function(evt) {
        var rect = canvas.getBoundingClientRect()
        if (reboundView.mouseDown==1){
            reboundView.mouseDown = 2;
            reboundView.mouse_x = evt.clientX-rect.left;
            reboundView.mouse_y = evt.clientY-rect.top;
            return;
        }else if (reboundView.mouseDown==2){
            var width = rect.width;
            var height = rect.height;
            var dx = 3.*(evt.clientX-rect.left-reboundView.mouse_x)/width;
            var dy = 3.*(evt.clientY-rect.top-reboundView.mouse_y)/height;
            reboundView.mouse_x = evt.clientX-rect.left;
            reboundView.mouse_y = evt.clientY-rect.top;
            if (evt.shiftKey){
                reboundView.scale *= (1.+dx+dy);
            }else{
                var inv = conjugate(reboundView.view);
                var up = [0.,1.,0.];
                var right = [1.,0.,0.];
                var inv_up = [0.,0.,0.];
                var inv_right = [0.,0.,0.];
                multvec(inv, right, inv_right);
                multvec(inv, up, inv_up);
                
                var sin_dy = Math.sin(dy);
                var rot_dy = {x:inv_right[0]*sin_dy, y:inv_right[1]*sin_dy, z:inv_right[2]*sin_dy, w:Math.cos(dy)};
                reboundView.view = mult(reboundView.view, normalize(rot_dy));
                
                var sin_dx = Math.sin(dx);
                var rot_dx = {x:inv_up[0]*sin_dx, y:inv_up[1]*sin_dx, z:inv_up[2]*sin_dx, w:Math.cos(dx)};
                reboundView.view = normalize(mult(reboundView.view, normalize(rot_dx)));
            }

            drawGL(reboundView);
        }


        }, false);

    reboundView.gl = canvas.getContext("webgl")||canvas.getContext("experimental-webgl");
    if (!reboundView.gl) {
      alert("Unable to initialize WebGL. Your browser may not support it.");
      return;
    }
    var gl = reboundView.gl
    gl.enable(gl.BLEND);
    gl.blendFunc(gl.SRC_ALPHA, gl.ONE_MINUS_SRC_ALPHA);
    
    reboundView.orbit_shader_program = createProgramFromScripts(gl,"orbit_shader-vs","orbit_shader-fs");
    reboundView.point_shader_program = createProgramFromScripts(gl,"point_shader-vs","point_shader-fs");
   
    var lintwopi = new Float32Array(500);
    for(i=0;i<500;i++){
        lintwopi[i] = 2.*Math.PI/500.*i;
    }
    reboundView.orbit_lintwopi_buffer = gl.createBuffer();
    gl.bindBuffer(gl.ARRAY_BUFFER, reboundView.orbit_lintwopi_buffer);
    gl.bufferData(gl.ARRAY_BUFFER, 4*500, gl.STATIC_DRAW);
    gl.bufferSubData(gl.ARRAY_BUFFER, 0, lintwopi)
    reboundView.orbit_shader_mvp_location = gl.getUniformLocation(reboundView.orbit_shader_program,"mvp");
    reboundView.orbit_shader_focus_location = gl.getUniformLocation(reboundView.orbit_shader_program,"focus");
    reboundView.orbit_shader_aef_location = gl.getUniformLocation(reboundView.orbit_shader_program,"aef");
    reboundView.orbit_shader_omegaOmegainc_location = gl.getUniformLocation(reboundView.orbit_shader_program,"omegaOmegainc");
    
    reboundView.particle_data_buffer = gl.createBuffer();
    gl.useProgram(reboundView.point_shader_program);
    reboundView.point_shader_mvp_location = gl.getUniformLocation(reboundView.point_shader_program,"mvp");
    
    updateRenderData(reboundView);
    gl.clearColor(0.0, 0.0, 0.0, 1.0);
    gl.clear(gl.COLOR_BUFFER_BIT);
    drawGL(reboundView);
}
function updateRenderData(reboundView){
    var overlay = document.getElementById("reboundoverlay-"+reboundView.cid);
    overlay.innerHTML = reboundView.model.get("overlay");
    var previousN = reboundView.N;
    reboundView.N = reboundView.model.get("N");
    reboundView.t = reboundView.model.get("t");
    reboundView.particle_data = reboundView.model.get('particle_data');
    if (reboundView.orbits){
        reboundView.orbit_data = reboundView.model.get('orbit_data');
    }
    var gl = reboundView.gl
    if (reboundView.N>0){
        gl.bindBuffer(gl.ARRAY_BUFFER, reboundView.particle_data_buffer);
        gl.bufferData(gl.ARRAY_BUFFER, reboundView.N*7*4, gl.DYNAMIC_DRAW);
        gl.bufferSubData(gl.ARRAY_BUFFER, 0, reboundView.particle_data)
    }
}
function drawGL(reboundView) {
    if (!reboundView.gl){
        return;
    }
    // Cleanup
    var gl = reboundView.gl
    gl.clearColor(0.0, 0.0, 0.0, 1.0);
    gl.clear(gl.COLOR_BUFFER_BIT);
    
    // Draw
    gl.useProgram(reboundView.point_shader_program);
    gl.bindBuffer(gl.ARRAY_BUFFER, reboundView.particle_data_buffer);
    var pvp = gl.getAttribLocation(reboundView.point_shader_program,"vp");
    gl.enableVertexAttribArray(pvp);
    gl.vertexAttribPointer(pvp, 3, gl.FLOAT, 0, 4*7,0); // 4 = size of float
    var projection = [0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.];
    if (reboundView.ratio>=1.){
        matortho(projection, 
                -1.6*reboundView.scale, 1.6*reboundView.scale,
                -1.6/reboundView.ratio*reboundView.scale, 1.6/reboundView.ratio*reboundView.scale,
                -2.5*reboundView.scale, 2.5*reboundView.scale);
    }else{
        matortho(projection, 
                -1.6*reboundView.ratio*reboundView.scale, 1.6*reboundView.ratio*reboundView.scale,
                -1.6*reboundView.scale, 1.6*reboundView.scale,
                -2.5*reboundView.scale, 2.5*reboundView.scale);
    }
    var view = [0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.];
    quat2mat(reboundView.view,view);
    var mvp = [0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.];
    matmult(projection,view,mvp);
    gl.uniformMatrix4fv(reboundView.point_shader_mvp_location,false,mattransp(mvp));
    gl.drawArrays(gl.POINTS,0,reboundView.N);
   
    if (reboundView.orbits){
        gl.useProgram(reboundView.orbit_shader_program);
        gl.bindBuffer(gl.ARRAY_BUFFER, reboundView.orbit_lintwopi_buffer);
        var ltp = gl.getAttribLocation(reboundView.orbit_shader_program,"lintwopi");
        gl.enableVertexAttribArray(ltp);
        gl.vertexAttribPointer(ltp, 1, gl.FLOAT, 0, 0,0); // 4 = size of float
        gl.uniformMatrix4fv(reboundView.orbit_shader_mvp_location,false,mattransp(mvp));

        // Need to do this one by one
        // because WebGL is not supporting
        // instancing:
        for(i=0;i<reboundView.N-1;i++){
            var focus = new Float32Array(reboundView.orbit_data.buffer,4*9*i,3);
            gl.uniform3fv(reboundView.orbit_shader_focus_location,focus);
            var aef = new Float32Array(reboundView.orbit_data.buffer,4*(9*i+3),3);
            gl.uniform3fv(reboundView.orbit_shader_aef_location,aef);
            var omegaOmegainc = new Float32Array(reboundView.orbit_data.buffer,4*(9*i+6),3);
            gl.uniform3fv(reboundView.orbit_shader_omegaOmegainc_location,omegaOmegainc);

            gl.drawArrays(gl.LINE_STRIP,0,500);
        }
    }
}
require.undef('rebound');
    define('rebound', ["@jupyter-widgets/base"], function(widgets) {
    var ReboundView = widgets.DOMWidgetView.extend({
        render: function() {
            this.el.innerHTML = '<span style="display: inline-block; position: relative;" width="'+this.model.get("width")+'" height="'+this.model.get("height")+'"><canvas style="border: none;" id="reboundcanvas-'+this.cid+'" width="'+this.model.get("width")+'" height="'+this.model.get("height")+'"></canvas><span style="position: absolute; color: #FFF; pointer-events:none;  bottom:5px; right:0px; padding-right:5px; font-family: monospace;" id="reboundoverlay-'+this.cid+'">REBOUND</span></span>';
            this.model.on('change:t', this.trigger_refresh, this);
            this.model.on('change:count', this.trigger_refresh, this);
            this.model.on('change:screenshotcount', this.take_screenshot, this);
            this.startCount = 0;
            this.gl = null;
            // Only copy those once
            this.scale = this.model.get("scale");
            this.width = this.model.get("width");
            this.height = this.model.get("height");
            this.orbits = this.model.get("orbits");
            this.orientation = this.model.get("orientation");
            startGL(this);
        },
        take_screenshot: function() {
            drawGL(this);
            var canvas = document.getElementById("reboundcanvas-"+this.cid);
            var img = canvas.toDataURL("image/png");
            this.model.set("screenshot",img, {updated_view: this});
            this.touch();
        },
        trigger_refresh: function() {
            updateRenderData(this);
            drawGL(this);
        },
    });
    return {
        ReboundView: ReboundView
    };
});
      
</script>
"""

import ipywidgets
ipywidgets_major_version = int((ipywidgets.__version__).split(".")[0])


if ipywidgets_major_version<7:
    js_code = js_code.replace("@jupyter-widgets/base", "jupyter-js-widgets")
    js_code = js_code.replace(".cid", ".id")

from ipywidgets import DOMWidget
import traitlets
import math
import base64
import sys
from ctypes import c_float, byref, create_string_buffer, c_int, c_char, pointer
from . import clibrebound
def savescreenshot(change):
    if len(change["new"]) and change["type"] =="change":
        w = change["owner"]
        bd = base64.b64decode(change["new"].split(",")[-1])
        if sys.version_info[0] < 3:
            with open(w.screenshotprefix+"%05d.png"%w.screenshotcountall, 'w') as f:
                f.write(bd)
        else:
            with open(w.screenshotprefix+"%05d.png"%w.screenshotcountall, 'bw') as f:
                f.write(bd)
        w.screenshotcountall += 1
        if len(w.times)>w.screenshotcount:
            nexttime = w.times[w.screenshotcount]
            if w.archive:
                sim = w.archive.getSimulation(w.times[w.screenshotcount],mode=w.mode)
                w.refresh(pointer(sim))
            else:
                w.simp.contents.integrate(w.times[w.screenshotcount])
            w.screenshotcount += 1
        else:
            w.unobserve(savescreenshot)
            w.times = None
            w.screenshotprefix = None

class Widget(DOMWidget):
    _view_name = traitlets.Unicode('ReboundView').tag(sync=True)
    _view_module = traitlets.Unicode('rebound').tag(sync=True)
    count = traitlets.Int(0).tag(sync=True)
    screenshotcount = traitlets.Int(0).tag(sync=True)
    t = traitlets.Float().tag(sync=True)
    N = traitlets.Int().tag(sync=True)
    overlay = traitlets.Unicode('REB WIdund').tag(sync=True)
    width = traitlets.Float().tag(sync=True)
    height = traitlets.Float().tag(sync=True)
    scale = traitlets.Float().tag(sync=True)
    particle_data = traitlets.CBytes(allow_none=True).tag(sync=True)
    orbit_data = traitlets.CBytes(allow_none=True).tag(sync=True)
    orientation = traitlets.Tuple().tag(sync=True)
    orbits = traitlets.Int().tag(sync=True)
    screenshot = traitlets.Unicode().tag(sync=True)
    def __init__(self,simulation,size=(200,200),orientation=(0.,0.,0.,1.),scale=None,autorefresh=True,orbits=True, overlay=True):
        """ 
        Initializes a Widget.

        Widgets provide real-time 3D interactive visualizations for REBOUND simulations 
        within Jupyter Notebooks. To use widgets, the ipywidgets package needs to be installed
        and enabled in your Jupyter notebook server. 

        Parameters
        ----------
        size : (int, int), optional
            Specify the size of the widget in pixels. The default is 200 times 200 pixels.
        orientation : (float, float, float, float), optional
            Specify the initial orientation of the view. The four floats correspond to the
            x, y, z, and w components of a quaternion. The quaternion will be normalized.
        scale : float, optional
            Set the initial scale of the view. If not set, the widget will determine the 
            scale automatically based on current particle positions.
        autorefresh : bool, optional
            The default value if True. The view is updated whenever a particle is added,
            removed and every 100th of a second while a simulation is running. If set 
            to False, then the user needs to manually call the refresh() function on the 
            widget. This might be useful if performance is an issue.
        orbits : bool, optional
            The default value for this is True and the widget will draw the instantaneous 
            orbits of the particles. For simulations in which particles are not on
            Keplerian orbits, the orbits shown will not be accurate. 
        overlay : string, optional
            Change the default text overlay. Set to None to hide all text.
        """
        self.screenshotcountall = 0
        self.width, self.height  = size
        self.t, self.N = simulation.t, simulation.N
        self.orientation = orientation
        self.autorefresh = autorefresh
        self.orbits = orbits
        self.useroverlay = overlay
        self.simp = pointer(simulation)
        clibrebound.reb_display_copy_data.restype = c_int
        if scale is None:
            self.scale = simulation.display_data.contents.scale
        else:
            self.scale = scale
        self.count += 1

        super(Widget, self).__init__()

    def refresh(self, simp=None, isauto=0):
        """ 
        Manually refreshes a widget.
        
        Note that this function can also be called using the wrapper function of
        the Simulation object: sim.refreshWidgets(). 
        """

        if simp==None:
            simp = self.simp
        if self.autorefresh==0 and isauto==1:
            return
        sim = simp.contents
        size_changed = clibrebound.reb_display_copy_data(simp)
        clibrebound.reb_display_prepare_data(simp,c_int(self.orbits))
        if sim.N>0:
            self.particle_data = (c_char * (4*7*sim.N)).from_address(sim.display_data.contents.particle_data).raw
            if self.orbits:
                self.orbit_data = (c_char * (4*9*(sim.N-1))).from_address(sim.display_data.contents.orbit_data).raw
        if size_changed:
            #TODO: Implement better GPU size change
            pass
        if self.useroverlay==True:
            self.overlay = "REBOUND (%s), N=%d, t=%g"%(sim.integrator,sim.N,sim.t)
        elif self.useroverlay==None or self.useroverlay==False:
            self.overlay = ""
        else:
            self.overlay = self.useroverlay + ", N=%d, t=%g"%(sim.N,sim.t)
        self.N = sim.N
        self.t = sim.t
        self.count += 1

    def takeScreenshot(self, times=None, prefix="./screenshot", resetCounter=False, archive=None,mode="snapshot"):
        """
        Take one or more screenshots of the widget and save the images to a file. 
        The images can be used to create a video.

        This function cannot be called multiple times within one cell.

        Note: this is a new feature and might not work on all systems.
        It was tested on python 2.7.10 and 3.5.2 on MacOSX.
        
        Parameters
        ----------
        times : (float, list), optional
            If this argument is not given a screenshot of the widget will be made 
            as it is (without integrating the simulation). If a float is given, then the
            simulation will be integrated to that time and then a screenshot will 
            be taken. If a list of floats is given, the simulation will be integrated
            to each time specified in the array. A separate screenshot for 
            each time will be saved.
        prefix : (str), optional
            This string will be part of the output filename for each image.
            Follow by a five digit integer and the suffix .png. By default the
            prefix is './screenshot' which outputs images in the current
            directory with the filnames screenshot00000.png, screenshot00001.png...
            Note that the prefix can include a directory.
        resetCounter : (bool), optional
            Resets the output counter to 0. 
        archive : (rebound.SimulationArchive), optional
            Use a REBOUND SimulationArchive. Thus, instead of integratating the 
            Simulation from the current time, it will use the SimulationArchive
            to load a snapshot. See examples for usage.
        mode : (string), optional
            Mode to use when querying the SimulationArchive. See SimulationArchive
            documentation for details. By default the value is "snapshot".

        Examples
        --------

        First, create a simulation and widget. All of the following can go in 
        one cell.

        >>> sim = rebound.Simulation()
        >>> sim.add(m=1.)
        >>> sim.add(m=1.e-3,x=1.,vy=1.)
        >>> w = sim.getWidget()
        >>> w

        The widget should show up. To take a screenshot, simply call 

        >>> w.takeScreenshot()

        A new file with the name screenshot00000.png will appear in the 
        current directory. 

        Note that the takeScreenshot command needs to be in a separate cell, 
        i.e. after you see the widget. 

        You can pass an array of times to the function. This allows you to 
        take multiple screenshots, for example to create a movie, 

        >>> times = [0,10,100]
        >>> w.takeScreenshot(times)

        """
        self.archive = archive
        if resetCounter:
            self.screenshotcountall = 0
        self.screenshotprefix = prefix
        self.screenshotcount = 0
        self.overlay = "REBOUND"
        self.screenshot = "" 
        if archive is None:
            if times is None:
                times = self.simp.contents.t
            try:
                # List
                len(times)
            except:
                # Float:
                times = [times]
            self.times = times
            self.observe(savescreenshot,names="screenshot")
            self.simp.contents.integrate(times[0])
            self.screenshotcount += 1 # triggers first screenshot
        else:
            if times is None:
                raise ValueError("Need times argument for archive mode.")
            try:
                len(times)
            except:
                raise ValueError("Need a list of times for archive mode.")
            self.times = times
            self.mode = mode
            self.observe(savescreenshot,names="screenshot")
            sim = archive.getSimulation(times[0],mode=mode)
            self.refresh(pointer(sim))
            self.screenshotcount += 1 # triggers first screenshot



            



    @staticmethod
    def getClientCode():
        return shader_code + js_code
