import numpy as np

X_MCP = 145
Y_MCP = 55
#ns limit need divided by 2
A_X_LEFT_NS  = -50.44742
A_X_RIGHT_NS = 54.91726

A_Y_LEFT_NS  = -34.71570
A_Y_RIGHT_NS = 23.59747

B_X_LEFT_NS  = -51.18360
B_X_RIGHT_NS = 53.40479

def init_param(SN_ID):
    #SN_ID 1 orbit 2 ground
    global B_Y_LEFT_NS, B_Y_RIGHT_NS, A_X_MAX_NS, A_Y_MAX_NS, A_NS_XLIM, A_NS_YLIM, A_X_CENTER_NS, A_Y_CENTER_NS, A_X_COF, A_Y_COF, B_X_MAX_NS, B_Y_MAX_NS, B_NS_XLIM, B_NS_YLIM, B_X_CENTER_NS, B_Y_CENTER_NS, B_X_COF, B_Y_COF
    if SN_ID == 1:
        #B_Y_LEFT_NS  = -16.68598
        #B_Y_RIGHT_NS = 41.62719
        B_Y_LEFT_NS  = -34.71570
        B_Y_RIGHT_NS = 23.59747

    if SN_ID == 2:
        B_Y_LEFT_NS  = -34.71570
        B_Y_RIGHT_NS = 23.59747

    A_X_MAX_NS = (abs(A_X_RIGHT_NS) + abs(A_X_LEFT_NS));
    A_Y_MAX_NS = (abs(A_Y_RIGHT_NS) + abs(A_Y_LEFT_NS));
    #need divided by 2 and yes
    A_NS_XLIM = A_X_MAX_NS/2;
    A_NS_YLIM = A_Y_MAX_NS/2;
    A_X_CENTER_NS = (A_X_LEFT_NS + A_X_RIGHT_NS);
    A_Y_CENTER_NS = (A_Y_LEFT_NS + A_Y_RIGHT_NS);
    A_X_COF = X_MCP/A_X_MAX_NS;
    A_Y_COF = Y_MCP/A_Y_MAX_NS;

    B_X_MAX_NS = abs(B_X_RIGHT_NS) + abs(B_X_LEFT_NS);
    B_Y_MAX_NS = abs(B_Y_RIGHT_NS) + abs(B_Y_LEFT_NS);
    B_NS_XLIM = B_X_MAX_NS/2;
    B_NS_YLIM = B_Y_MAX_NS/2;
    #should divided by 2 but not yet
    B_X_CENTER_NS = (B_X_LEFT_NS + B_X_RIGHT_NS);
    B_Y_CENTER_NS = (B_Y_LEFT_NS + B_Y_RIGHT_NS);
    B_X_COF = X_MCP/B_X_MAX_NS;
    B_Y_COF = Y_MCP/B_Y_MAX_NS;



if __name__ == '__main__':
    print('test')