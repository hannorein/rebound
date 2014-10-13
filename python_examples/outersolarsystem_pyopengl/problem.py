import array
import numpy as np
from OpenGL.GL import *
from OpenGL.GLUT import *
from OpenGL.GLU import *
from OpenGL.constants import GLfloat_3,GLfloat_4


# Import the rebound module
import sys
sys.path.append('../')
import rebound
from rebound import Particle

################################################################################
## Global variables
################################################################################
visual_simulation = True
# timestep counter
steps = 0
display_wire = True
# Window
g_fViewDistance = 100.
g_Width = 700
g_Height = 700
g_nearPlane = 10
g_farPlane = 1000
# Movement
action = ""
xStart = yStart = 0.
zoom = 100.
xRotate = 0.
yRotate = 0.
zRotate = 0.
xTrans = 0.
yTrans = 0.
################################################################################

def mouse(button, state, x, y):
    global action, xStart, yStart
    if (button==GLUT_LEFT_BUTTON):
        if (glutGetModifiers() == GLUT_ACTIVE_SHIFT):
            action = "MOVE_EYE_2"
        else:
            action = "MOVE_EYE"
    elif (button==GLUT_MIDDLE_BUTTON):
        action = "ZOOM"
    elif (button==GLUT_RIGHT_BUTTON):
        action = "TRANS"
    xStart = x
    yStart = y


def motion(x, y):
    global zoom, xStart, yStart, xRotate, yRotate, zRotate, xTrans, yTrans
    if (action=="MOVE_EYE"):
        xRotate += x - xStart
        yRotate -= y - yStart
    elif (action=="MOVE_EYE_2"):
        zRotate += y - yStart
    elif (action=="TRANS"):
        xTrans += x - xStart
        yTrans += y - yStart
    elif (action=="ZOOM"):
        zoom -= y - yStart
        if zoom > 150.:
            zoom = 150.
        elif zoom < 1.1:
            zoom = 1.1
    else:
        print("unknown action\n", action)
    xStart = x
    yStart = y
    glutPostRedisplay()


def resetView():
    global zoom, xRotate, yRotate, zRotate, xTrans, yTrans
    zoom = 65.
    xRotate = 0.
    yRotate = 0.
    zRotate = 0.
    xTrans = 0.
    yTrans = 0.
    glutPostRedisplay()


def reshape(width, height):
    global g_Width, g_Height, g_nearPlane, g_farPlane
    g_Width = width
    g_Height = height
    glViewport(0, 0, g_Width, g_Height)

    _top    =  1.0
    _bottom = -1.0
    _left   = -g_Width/g_Height
    _right  = -_left
    _zNear  = -10.0
    _zFar   = 10.0

    glMatrixMode(GL_PROJECTION)
    glLoadIdentity()
    #  Orthographic Projection does not apply a perspective distortion and the perspective one does.
    #  Perspective Projection is the natural one as we Humans see in a perspective Way
    #  (things farther away from us appear smaller)
    #glOrtho(_left,_right,_bottom,_top,_zNear,_zFar)
    gluPerspective(zoom, float(g_Width)/float(g_Height), g_nearPlane, g_farPlane)
    glMatrixMode(GL_MODELVIEW)

def keyboard(key, x, y):
    if(key=='r'): resetView()
    if(key=='q'): exit(0)
    glutPostRedisplay()


def scenemodel():
    global particles, N, N_active

    if not display_wire:
        glEnable(GL_BLEND)
        glDepthMask(GL_FALSE)
        glDisable(GL_DEPTH_TEST)
        glDisable(GL_LIGHTING)
        glDisable(GL_LIGHT0)
    glEnable(GL_POINT_SMOOTH)

    vertices = array.array('f', [])
    for n in xrange(N):
        vertices += array.array('f', [particles[n].x, particles[n].y, particles[n].z])

    glVertexPointer( 3, GL_FLOAT, 0, vertices.tostring( ) )

    # Draw points
    glEnableClientState(GL_VERTEX_ARRAY)
    glPointSize(3.)
    glColor4f(1.0,1.0,1.0,0.5)
    glDrawArrays(GL_POINTS, N_active, N-N_active)
    glColor4f(1.0,1.0,0.0,0.9)
    glPointSize(5.)
    glDrawArrays(GL_POINTS, 0, N_active)
    glDisableClientState(GL_VERTEX_ARRAY)

    # Draw orbits
    if display_wire:
        radius = 0
        com = particles[0]
        for i in range(1, N):
            p = particles[i]
            if (N_active>0):
                #// Different colors for active/test particles
                if (i>=N_active):
                    glColor4f(0.9,1.0,0.9,0.9)
                else:
                    glColor4f(1.0,0.9,0.0,0.9)
            else:
                #// Alternating colors
                if (i%2 == 1):
                    glColor4f(0.0,1.0,0.0,0.9)
                else:
                    glColor4f(0.0,0.0,1.0,0.9)
            o = rebound.tools_p2orbit(p,com)
            glPushMatrix()

            glTranslatef(com.x,com.y,com.z)
            DEG2RAD = np.pi/180
            glRotatef(o.Omega/DEG2RAD,0,0,1)
            glRotatef(o.inc/DEG2RAD,1,0,0)
            glRotatef(o.omega/DEG2RAD,0,0,1)

            glBegin(GL_LINE_LOOP)
            for trueAnom in np.arange(0, 2*np.pi, np.pi/100.):
                #//convert degrees into radians
                radius = o.a * (1. - o.e*o.e) / (1. + o.e*np.cos(trueAnom))
                glVertex3f(radius*np.cos(trueAnom),radius*np.sin(trueAnom),0)
            glEnd()
            glPopMatrix()
            com = rebound.tools_get_center_of_mass(p,com)


    #glMatrixMode(GL_PROJECTION)
    boxsize_x = 100
    boxsize_y = 100
    boxsize_z = 50
    glColor4f(1.0,0.0,0.0,0.4)
    glScalef(boxsize_x,boxsize_y,boxsize_z)
    glutWireCube(1)
    glScalef(1./boxsize_x,1./boxsize_y,1./boxsize_z)
    #glMatrixMode(GL_MODELVIEW) # the "default" and safest mode to leave OpenGL in

def polarView():
    glTranslatef( yTrans/100., 0.0, 0.0 )
    glTranslatef(  0.0, -xTrans/100., 0.0)
    glRotatef( -zRotate, 0.0, 0.0, 1.0)
    glRotatef( -xRotate, 1.0, 0.0, 0.0)
    glRotatef( -yRotate, .0, 1.0, 0.0)


def display():
    global g_nearPlane, g_farPlane, g_Height, g_Width, g_fViewDistance
    # Add a projection matrix (perspective)
    glMatrixMode(GL_PROJECTION)
    glLoadIdentity() # Reset the projection matrix, or bad things happen after multiple calls to below functions!
    #glOrtho( ... ) // to use 2D rendering
    #gluPerspective( ... ) // to use (easy) 3D rendering
    #glFrustrum( ... ) // to use (harder/less intuitive?) 3D rendering
    gluPerspective(zoom, float(g_Width)/float(g_Height), g_nearPlane, g_farPlane)
    #   - field of view angle, in degrees, in the y direction
    #   - aspect ratio that determines the field of view in the x direction
    #   - distance from the viewer to the near clipping plane
    #   - distance from the viewer to the far clipping plane
    #       * The greater the ratio of zFar to zNear is, the less effective the depth buffer will be at distinguishing between surfaces that are near each other
    glMatrixMode(GL_MODELVIEW) # the "default" and safest mode to leave OpenGL in

    # Clear frame buffer and depth buffer
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
    glLoadIdentity() # Reset the projection matrix
    # Set camera / geometry position (looking down -Z axis)
    gluLookAt(0, 0, -g_fViewDistance, 0, 0, 0, -.1, 0, 0)   # Eye X Y Z, reference point X Y Z, up vector X Y Z

    # Render the scene
    polarView()
    scenemodel()
    # Make sure changes appear onscreen
    glutSwapBuffers()


def visualize(iterate_func):
    global visual_simulation
    visual_simulation = True
    ############################################################################
    ## OpenGL
    ############################################################################
    glutInit()
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH )
    glutInitWindowSize (g_Width,g_Height)
    glutInitWindowPosition (0 + 4, g_Height / 4)
    glutCreateWindow("rebound")
    # Register callbacks
    glutReshapeFunc(reshape)
    glutMouseFunc(mouse)
    glutMotionFunc(motion)
    glutDisplayFunc(display)
    glutIdleFunc(iterate_func)
    glutKeyboardFunc(keyboard)
    # Enter glut run loop and never come back.
    glutMainLoop()
    ############################################################################


def iterate():
    global steps, visual_simulation
    rebound.step()
    steps += 1
    # Print particle positions every 100 timesteps
    if steps%1000==0:
        print rebound.get_t()
        for i in range(rebound.get_N()):
            #     time             particle id   x               y               z
            print "\t", i,            particles[i].x, particles[i].y, particles[i].z

    if visual_simulation:
        display()


if __name__ == '__main__':
    # Set variables (defaults are G=1, t=0, dt=0.01)
    k = 0.01720209895       # Gaussian constant
    rebound.set_G(k*k)      # Gravitational constant

    # Setup particles (data taken from NASA Horizons)
    # This could also be easily read in from a file.
    rebound.particle_add( Particle( m=1.00000597682, x=-4.06428567034226e-3, y=-6.08813756435987e-3, z=-1.66162304225834e-6, vx=+6.69048890636161e-6, vy=-6.33922479583593e-6, vz=-3.13202145590767e-9) )  # Sun
    rebound.particle_add( Particle( m=1./1047.355,   x=+3.40546614227466e+0, y=+3.62978190075864e+0, z=+3.42386261766577e-2, vx=-5.59797969310664e-3, vy=+5.51815399480116e-3, vz=-2.66711392865591e-6) )  # Jupiter
    rebound.particle_add( Particle( m=1./3501.6,     x=+6.60801554403466e+0, y=+6.38084674585064e+0, z=-1.36145963724542e-1, vx=-4.17354020307064e-3, vy=+3.99723751748116e-3, vz=+1.67206320571441e-5) )  # Saturn
    rebound.particle_add( Particle( m=1./22869.,     x=+1.11636331405597e+1, y=+1.60373479057256e+1, z=+3.61783279369958e-1, vx=-3.25884806151064e-3, vy=+2.06438412905916e-3, vz=-2.17699042180559e-5) )  # Uranus
    rebound.particle_add( Particle( m=1./19314.,     x=-3.01777243405203e+1, y=+1.91155314998064e+0, z=-1.53887595621042e-1, vx=-2.17471785045538e-4, vy=-3.11361111025884e-3, vz=+3.58344705491441e-5) )  # Neptune
    rebound.particle_add( Particle( m=0,             x=-2.13858977531573e+1, y=+3.20719104739886e+1, z=+2.49245689556096e+0, vx=-1.76936577252484e-3, vy=-2.06720938381724e-3, vz=+6.58091931493844e-4) )  # Pluto
    N = rebound.get_N()
    N_active = 5

    # Set the center of momentum to be at the origin
    rebound.move_to_center_of_momentum()

    # Get the particle data
    # Note: this is a pointer and will automatically update as the simulation progresses
    particles = rebound.particles_get()
    orbit = rebound.tools_p2orbit(particles[1], particles[0])
    if visual_simulation:
        visualize(iterate)
    else:
        # Integrate until t=1e6 (unit of time in this example is days)
        while rebound.get_t()<1e6:
            iterate()


