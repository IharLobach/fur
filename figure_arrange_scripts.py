import os
from PIL import Image,ImageDraw,ImageFont
from IPython.display import display


def RowOfFigures(im_paths,output_path,hor_separ=0.05,v_space=100,fontsize=200,label_pos_x=0,label_pos_y=0,display_inline=1):
    """Creates an image with a row of figures labeled 'a','b','c', etc. in the top-left corner
    Parameters
    ----------
    im_paths
        list of input images' paths
    output_path
        path where to save the output image
    hor_separ
        horizontal separation in units of sum of widths of all figures
    v_space : int
        pixels, additional vertical white space on top of the figure (helps to fit the labels 'a','b', etc.)
    fontsize
        labels' fontsize
    label_pos_x : int
        x position of a label
    label_pos_y :int
        y position of a label
    display_inline
        0 or 1, 0 means not displayed, 1 means displayed inline in the notebook
    """
    images = [Image.open(p) for p in im_paths]
    ws = []
    hs = []
    for i in images:
        ws.append(i.size[0])
        hs.append(i.size[1])
    for w in ws:
        if w!=ws[0]:
            print('Warning: Widths are not equal!')
            break
    for h in hs:
        if h!=hs[0]:
            print('Warning: Heights are not equal!')
            break
    wmax = max(ws)
    hmax = max(hs)
    wsum = sum(ws)
    gap = round(hor_separ*wsum)
    width = wsum +(len(ws)-1)*gap
    height = hmax+v_space
    result = Image.new("RGB", (width,height),color=(255,255,255,0))
    d = ImageDraw.Draw(result)
    # use a truetype font
    font = ImageFont.truetype(
        "/usr/share/fonts/truetype/msttcorefonts/arial.ttf", fontsize)
    x=0
    q=0
    labels=['a','b','c','d','e','f','g']
    for i in images:
        w,h=i.size
        result.paste(i,(x,height-h,x+w,height))
        d.text((x+label_pos_x,label_pos_y),labels[q],fill=(0,0,0) ,font=font)
        x+=w+gap
        q+=1
    result.save(output_path)
    print('Result saved to '+output_path)
    if display_inline==1:
        display(result)

def FourFigures2x2(im_paths,output_path,hor_separ=0,ver_separ=0,v_space=0,fontsize=20,label_pos_x=0,label_pos_y=0,display_inline=0):
    """Creates an image with a row of figures labeled 'a','b','c', etc. in the top-left corner
    Parameters
    ----------
    im_paths
        list of input images' paths
    output_path
        path where to save the output image
    hor_separ
        horizontal separation in units of sum of widths of all figures
    ver_separ : int
        pixels, vertical separation between the two rows of subfigures
    v_space : int
        pixels, additional vertical white space on top of the figure (helps to fit the labels 'a','b', etc.)
    fontsize
        labels' fontsize
    label_pos_x : int
        x position of a label
    label_pos_y :int
        y position of a label
    display_inline
        0 or 1, 0 means not displayed, 1 means displayed inline in the notebook
    """
    images = [Image.open(p) for p in im_paths]
    ws = []
    hs = []
    for i in images:
        ws.append(i.size[0])
        hs.append(i.size[1])
    print("In this script figures' widths must be equal!")
    for w in ws:
        if w!=ws[0]:
            print('Error: Widths are not equal!')
            break

    for h in hs:
        if h!=hs[0]:
            print("Warning: Heights are not equal! But it's okay.")
            break
    hmax1 = max(hs[0:2])
    hmax2 = max(hs[2:4])
    wsum = 2*ws[0]
    hgap = round(hor_separ*wsum)
    width = wsum +hgap
    height = hmax1+2*v_space+hmax2+ver_separ
    result = Image.new("RGB", (width,height),color=(255,255,255,0))
    d = ImageDraw.Draw(result)
    # use a truetype font
    font = ImageFont.truetype("arial.ttf", fontsize)
    labels1=['a','b']
    labels2=['c','d']
    x=0
    q=0
    for i in images[2:4]:
        w,h=i.size
        result.paste(i,(x,height-h,x+w,height))
        d.text((x+label_pos_x,label_pos_y+hmax1+v_space+ver_separ),labels2[q],fill=(0,0,0),font=font)
        x+=w+hgap
        q+=1
    x=0
    q=0
    for i in images[0:2]:
        w,h=i.size
        result.paste(i,(x,hmax1+v_space-h,x+w,hmax1+v_space))
        d.text((x+label_pos_x,label_pos_y),labels1[q],fill=(0,0,0),font=font)
        x+=w+hgap
        q+=1
    result.save(output_path)
    print('Result saved to '+output_path)
    if display_inline==1:
        display(result)

def as_si(x, ndp):
    """Returns number in scientific format for latex. For PyPlot when latex output is needed. 
    Latex distibtuion must be installed on the machine. Also set  plt.rcParams["text.usetex"] =True
    Parameters
    ----------
    x
        the input number
    ndp : int
        desired number of decimals
    """
    s = '{x:0.{ndp:d}e}'.format(x=x, ndp=ndp)
    m, e = s.split('e')
    return r'${m:s}\times 10^{{{e:d}}}$'.format(m=m, e=int(e))
