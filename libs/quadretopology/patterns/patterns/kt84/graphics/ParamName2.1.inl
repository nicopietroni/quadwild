// returns one value, the number of alpha bitplanes in the accumulation buffer.
ACCUM_ALPHA_BITS = GL_ACCUM_ALPHA_BITS,
// returns one value, the number of blue bitplanes in the accumulation buffer.
ACCUM_BLUE_BITS = GL_ACCUM_BLUE_BITS,
// returns four values: the red, green, blue, and alpha values used to clear the accumulation buffer. Integer values, if requested, are linearly mapped from the internal floating-point representation such that 1.0 returns the most positive representable integer value, and -1.0 returns the most negative representable integer value. The initial value is (0, 0, 0, 0). See glClearAccum.
ACCUM_CLEAR_VALUE = GL_ACCUM_CLEAR_VALUE,
// returns one value, the number of green bitplanes in the accumulation buffer.
ACCUM_GREEN_BITS = GL_ACCUM_GREEN_BITS,
// returns one value, the number of red bitplanes in the accumulation buffer.
ACCUM_RED_BITS = GL_ACCUM_RED_BITS,
// returns a single value indicating the active multitexture unit. The initial value is GL_TEXTURE0. See glActiveTexture.
ACTIVE_TEXTURE = GL_ACTIVE_TEXTURE,
// returns two values, the smallest and largest supported sizes for aliased points.
ALIASED_POINT_SIZE_RANGE = GL_ALIASED_POINT_SIZE_RANGE,
// returns two values, the smallest and largest supported widths for aliased lines.
ALIASED_LINE_WIDTH_RANGE = GL_ALIASED_LINE_WIDTH_RANGE,
// returns one value, the alpha bias factor used during pixel transfers. The initial value is 0. See glPixelTransfer.
ALPHA_BIAS = GL_ALPHA_BIAS,
// returns one value, the number of alpha bitplanes in each color buffer.
ALPHA_BITS = GL_ALPHA_BITS,
// returns one value, the alpha scale factor used during pixel transfers. The initial value is 1. See glPixelTransfer.
ALPHA_SCALE = GL_ALPHA_SCALE,
// returns a single boolean value indicating whether alpha testing of fragments is enabled. The initial value is GL_FALSE. See glAlphaFunc.
ALPHA_TEST = GL_ALPHA_TEST,
// returns one value, the symbolic name of the alpha test function. The initial value is GL_ALWAYS. See glAlphaFunc.
ALPHA_TEST_FUNC = GL_ALPHA_TEST_FUNC,
// returns one value, the reference value for the alpha test. The initial value is 0. See glAlphaFunc. An integer value, if requested, is linearly mapped from the internal floating-point representation such that 1.0 returns the most positive representable integer value, and -1.0 returns the most negative representable integer value.
ALPHA_TEST_REF = GL_ALPHA_TEST_REF,
// returns a single value, the name of the buffer object currently bound to the target GL_ARRAY_BUFFER. If no buffer object is bound to this target, 0 is returned. The initial value is 0. See glBindBuffer.
ARRAY_BUFFER_BINDING = GL_ARRAY_BUFFER_BINDING,
// returns one value, the depth of the attribute stack. If the stack is empty, 0 is returned. The initial value is 0. See glPushAttrib.
ATTRIB_STACK_DEPTH = GL_ATTRIB_STACK_DEPTH,
// returns a single boolean value indicating whether 2D map evaluation automatically generates surface normals. The initial value is GL_FALSE. See glMap2.
AUTO_NORMAL = GL_AUTO_NORMAL,
// returns one value, the number of auxiliary color buffers available.
AUX_BUFFERS = GL_AUX_BUFFERS,
// returns a single boolean value indicating whether blending is enabled. The initial value is GL_FALSE. See glBlendFunc.
BLEND = GL_BLEND,
// returns four values, the red, green, blue, and alpha values which are the components of the blend color. See glBlendColor.
BLEND_COLOR = GL_BLEND_COLOR,
// returns one value, the symbolic constant identifying the alpha destination blend function. The initial value is GL_ZERO. See glBlendFunc and glBlendFuncSeparate.
BLEND_DST_ALPHA = GL_BLEND_DST_ALPHA,
// returns one value, the symbolic constant identifying the RGB destination blend function. The initial value is GL_ZERO. See glBlendFunc and glBlendFuncSeparate.
BLEND_DST_RGB = GL_BLEND_DST_RGB,
// returns one value, a symbolic constant indicating whether the RGB blend equation is GL_FUNC_ADD, GL_FUNC_SUBTRACT, GL_FUNC_REVERSE_SUBTRACT, GL_MIN or GL_MAX. See glBlendEquationSeparate.
BLEND_EQUATION_RGB = GL_BLEND_EQUATION_RGB,
// returns one value, a symbolic constant indicating whether the Alpha blend equation is GL_FUNC_ADD, GL_FUNC_SUBTRACT, GL_FUNC_REVERSE_SUBTRACT, GL_MIN or GL_MAX. See glBlendEquationSeparate.
BLEND_EQUATION_ALPHA = GL_BLEND_EQUATION_ALPHA,
// returns one value, the symbolic constant identifying the alpha source blend function. The initial value is GL_ONE. See glBlendFunc and glBlendFuncSeparate.
BLEND_SRC_ALPHA = GL_BLEND_SRC_ALPHA,
// returns one value, the symbolic constant identifying the RGB source blend function. The initial value is GL_ONE. See glBlendFunc and glBlendFuncSeparate.
BLEND_SRC_RGB = GL_BLEND_SRC_RGB,
// returns one value, the blue bias factor used during pixel transfers. The initial value is 0. See glPixelTransfer.
BLUE_BIAS = GL_BLUE_BIAS,
// returns one value, the number of blue bitplanes in each color buffer.
BLUE_BITS = GL_BLUE_BITS,
// returns one value, the blue scale factor used during pixel transfers. The initial value is 1. See glPixelTransfer.
BLUE_SCALE = GL_BLUE_SCALE,
// returns a single integer value indicating the current client active multitexture unit. The initial value is GL_TEXTURE0. See glClientActiveTexture.
CLIENT_ACTIVE_TEXTURE = GL_CLIENT_ACTIVE_TEXTURE,
// returns one value indicating the depth of the attribute stack. The initial value is 0. See glPushClientAttrib.
CLIENT_ATTRIB_STACK_DEPTH = GL_CLIENT_ATTRIB_STACK_DEPTH,
// returns a single boolean value indicating whether the specified clipping plane is enabled. The initial value is GL_FALSE. See glClipPlane.
CLIP_PLANE0 = GL_CLIP_PLANE0,
CLIP_PLANE1 = GL_CLIP_PLANE1,
CLIP_PLANE2 = GL_CLIP_PLANE2,
CLIP_PLANE3 = GL_CLIP_PLANE3,
CLIP_PLANE4 = GL_CLIP_PLANE4,
CLIP_PLANE5 = GL_CLIP_PLANE5,
// returns a single boolean value indicating whether the color array is enabled. The initial value is GL_FALSE. See glColorPointer.
COLOR_ARRAY = GL_COLOR_ARRAY,
// returns a single value, the name of the buffer object associated with the color array. This buffer object would have been bound to the target GL_ARRAY_BUFFER at the time of the most recent call to glColorPointer. If no buffer object was bound to this target, 0 is returned. The initial value is 0. See glBindBuffer.
COLOR_ARRAY_BUFFER_BINDING = GL_COLOR_ARRAY_BUFFER_BINDING,
// returns one value, the number of components per color in the color array. The initial value is 4. See glColorPointer.
COLOR_ARRAY_SIZE = GL_COLOR_ARRAY_SIZE,
// returns one value, the byte offset between consecutive colors in the color array. The initial value is 0. See glColorPointer.
COLOR_ARRAY_STRIDE = GL_COLOR_ARRAY_STRIDE,
// returns one value, the data type of each component in the color array. The initial value is GL_FLOAT. See glColorPointer.
COLOR_ARRAY_TYPE = GL_COLOR_ARRAY_TYPE,
// returns four values: the red, green, blue, and alpha values used to clear the color buffers. Integer values, if requested, are linearly mapped from the internal floating-point representation such that 1.0 returns the most positive representable integer value, and -1.0 returns the most negative representable integer value. The initial value is (0, 0, 0, 0). See glClearColor.
COLOR_CLEAR_VALUE = GL_COLOR_CLEAR_VALUE,
// returns a single boolean value indicating whether a fragment's RGBA color values are merged into the framebuffer using a logical operation. The initial value is GL_FALSE. See glLogicOp.
COLOR_LOGIC_OP = GL_COLOR_LOGIC_OP,
// returns a single boolean value indicating whether one or more material parameters are tracking the current color. The initial value is GL_FALSE. See glColorMaterial.
COLOR_MATERIAL = GL_COLOR_MATERIAL,
// returns one value, a symbolic constant indicating which materials have a parameter that is tracking the current color. The initial value is GL_FRONT_AND_BACK. See glColorMaterial.
COLOR_MATERIAL_FACE = GL_COLOR_MATERIAL_FACE,
// returns one value, a symbolic constant indicating which material parameters are tracking the current color. The initial value is GL_AMBIENT_AND_DIFFUSE. See glColorMaterial.
COLOR_MATERIAL_PARAMETER = GL_COLOR_MATERIAL_PARAMETER,
// returns sixteen values: the color matrix on the top of the color matrix stack. Initially this matrix is the identity matrix. See glPushMatrix.
COLOR_MATRIX = GL_COLOR_MATRIX,
// returns one value, the maximum supported depth of the projection matrix stack. The value must be at least 2. See glPushMatrix.
COLOR_MATRIX_STACK_DEPTH = GL_COLOR_MATRIX_STACK_DEPTH,
// returns a single boolean value indicating whether primary and secondary color sum is enabled. See glSecondaryColor.
COLOR_SUM = GL_COLOR_SUM,
// returns a single boolean value indicating whether the color table lookup is enabled. See glColorTable.
COLOR_TABLE = GL_COLOR_TABLE,
// returns four boolean values: the red, green, blue, and alpha write enables for the color buffers. The initial value is (GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE). See glColorMask.
COLOR_WRITEMASK = GL_COLOR_WRITEMASK,
// returns a list of symbolic constants of length GL_NUM_COMPRESSED_TEXTURE_FORMATS indicating which compressed texture formats are available. See glCompressedTexImage2D.
COMPRESSED_TEXTURE_FORMATS = GL_COMPRESSED_TEXTURE_FORMATS,
// returns a single boolean value indicating whether 1D convolution is enabled. The initial value is GL_FALSE. See glConvolutionFilter1D.
CONVOLUTION_1D = GL_CONVOLUTION_1D,
// returns a single boolean value indicating whether 2D convolution is enabled. The initial value is GL_FALSE. See glConvolutionFilter2D.
CONVOLUTION_2D = GL_CONVOLUTION_2D,
// returns a single boolean value indicating whether polygon culling is enabled. The initial value is GL_FALSE. See glCullFace.
CULL_FACE = GL_CULL_FACE,
// returns one value, a symbolic constant indicating which polygon faces are to be culled. The initial value is GL_BACK. See glCullFace.
CULL_FACE_MODE = GL_CULL_FACE_MODE,
// returns four values: the red, green, blue, and alpha values of the current color. Integer values, if requested, are linearly mapped from the internal floating-point representation such that 1.0 returns the most positive representable integer value, and -1.0 returns the most negative representable integer value. The initial value is (1, 1, 1, 1). See glColor.
CURRENT_COLOR = GL_CURRENT_COLOR,
// returns one value, the current fog coordinate. The initial value is 0. See glFogCoord.
CURRENT_FOG_COORD = GL_CURRENT_FOG_COORD,
// returns one value, the current color index. The initial value is 1. See glIndex.
CURRENT_INDEX = GL_CURRENT_INDEX,
// returns three values: the x, y, and z values of the current normal. Integer values, if requested, are linearly mapped from the internal floating-point representation such that 1.0 returns the most positive representable integer value, and -1.0 returns the most negative representable integer value. The initial value is (0, 0, 1). See glNormal.
CURRENT_NORMAL = GL_CURRENT_NORMAL,
// returns one value, the name of the program object that is currently active, or 0 if no program object is active. See glUseProgram.
CURRENT_PROGRAM = GL_CURRENT_PROGRAM,
// returns four values: the red, green, blue, and alpha color values of the current raster position. Integer values, if requested, are linearly mapped from the internal floating-point representation such that 1.0 returns the most positive representable integer value, and -1.0 returns the most negative representable integer value. The initial value is (1, 1, 1, 1). See glRasterPos.
CURRENT_RASTER_COLOR = GL_CURRENT_RASTER_COLOR,
// returns one value, the distance from the eye to the current raster position. The initial value is 0. See glRasterPos.
CURRENT_RASTER_DISTANCE = GL_CURRENT_RASTER_DISTANCE,
// returns one value, the color index of the current raster position. The initial value is 1. See glRasterPos.
CURRENT_RASTER_INDEX = GL_CURRENT_RASTER_INDEX,
// returns four values: the x, y, z, and w components of the current raster position. x, y, and z are in window coordinates, and w is in clip coordinates. The initial value is (0, 0, 0, 1). See glRasterPos.
CURRENT_RASTER_POSITION = GL_CURRENT_RASTER_POSITION,
// returns a single boolean value indicating whether the current raster position is valid. The initial value is GL_TRUE. See glRasterPos.
CURRENT_RASTER_POSITION_VALID = GL_CURRENT_RASTER_POSITION_VALID,
// returns four values: the red, green, blue, and alpha secondary color values of the current raster position. Integer values, if requested, are linearly mapped from the internal floating-point representation such that 1.0 returns the most positive representable integer value, and -1.0 returns the most negative representable integer value. The initial value is (1, 1, 1, 1). See glRasterPos.
CURRENT_RASTER_SECONDARY_COLOR = GL_CURRENT_RASTER_SECONDARY_COLOR,
// returns four values: the s, t, r, and q texture coordinates of the current raster position. The initial value is (0, 0, 0, 1). See glRasterPos and glMultiTexCoord.
CURRENT_RASTER_TEXTURE_COORDS = GL_CURRENT_RASTER_TEXTURE_COORDS,
// returns four values: the red, green, blue, and alpha values of the current secondary color. Integer values, if requested, are linearly mapped from the internal floating-point representation such that 1.0 returns the most positive representable integer value, and -1.0 returns the most negative representable integer value. The initial value is (0, 0, 0, 0). See glSecondaryColor.
CURRENT_SECONDARY_COLOR = GL_CURRENT_SECONDARY_COLOR,
// returns four values: the s, t, r, and q current texture coordinates. The initial value is (0, 0, 0, 1). See glMultiTexCoord.
CURRENT_TEXTURE_COORDS = GL_CURRENT_TEXTURE_COORDS,
// returns one value, the depth bias factor used during pixel transfers. The initial value is 0. See glPixelTransfer.
DEPTH_BIAS = GL_DEPTH_BIAS,
// returns one value, the number of bitplanes in the depth buffer.
DEPTH_BITS = GL_DEPTH_BITS,
// returns one value, the value that is used to clear the depth buffer. Integer values, if requested, are linearly mapped from the internal floating-point representation such that 1.0 returns the most positive representable integer value, and -1.0 returns the most negative representable integer value. The initial value is 1. See glClearDepth.
DEPTH_CLEAR_VALUE = GL_DEPTH_CLEAR_VALUE,
// returns one value, the symbolic constant that indicates the depth comparison function. The initial value is GL_LESS. See glDepthFunc.
DEPTH_FUNC = GL_DEPTH_FUNC,
// returns two values: the near and far mapping limits for the depth buffer. Integer values, if requested, are linearly mapped from the internal floating-point representation such that 1.0 returns the most positive representable integer value, and -1.0 returns the most negative representable integer value. The initial value is (0, 1). See glDepthRange.
DEPTH_RANGE = GL_DEPTH_RANGE,
// returns one value, the depth scale factor used during pixel transfers. The initial value is 1. See glPixelTransfer.
DEPTH_SCALE = GL_DEPTH_SCALE,
// returns a single boolean value indicating whether depth testing of fragments is enabled. The initial value is GL_FALSE. See glDepthFunc and glDepthRange.
DEPTH_TEST = GL_DEPTH_TEST,
// returns a single boolean value indicating if the depth buffer is enabled for writing. The initial value is GL_TRUE. See glDepthMask.
DEPTH_WRITEMASK = GL_DEPTH_WRITEMASK,
// returns a single boolean value indicating whether dithering of fragment colors and indices is enabled. The initial value is GL_TRUE.
DITHER = GL_DITHER,
// returns a single boolean value indicating whether double buffering is supported.
DOUBLEBUFFER = GL_DOUBLEBUFFER,
// returns one value, a symbolic constant indicating which buffers are being drawn to. See glDrawBuffer. The initial value is GL_BACK if there are back buffers, otherwise it is GL_FRONT.
DRAW_BUFFER = GL_DRAW_BUFFER,
// returns one value, a symbolic constant indicating which buffers are being drawn to by the corresponding output color. See glDrawBuffers. The initial value of GL_DRAW_BUFFER0 is GL_BACK if there are back buffers, otherwise it is GL_FRONT. The initial values of draw buffers for all other output colors is GL_NONE.
DRAW_BUFFER0  = GL_DRAW_BUFFER0 ,
DRAW_BUFFER1  = GL_DRAW_BUFFER1 ,
DRAW_BUFFER2  = GL_DRAW_BUFFER2 ,
DRAW_BUFFER3  = GL_DRAW_BUFFER3 ,
DRAW_BUFFER4  = GL_DRAW_BUFFER4 ,
DRAW_BUFFER5  = GL_DRAW_BUFFER5 ,
DRAW_BUFFER6  = GL_DRAW_BUFFER6 ,
DRAW_BUFFER7  = GL_DRAW_BUFFER7 ,
DRAW_BUFFER8  = GL_DRAW_BUFFER8 ,
DRAW_BUFFER9  = GL_DRAW_BUFFER9 ,
DRAW_BUFFER10 = GL_DRAW_BUFFER10,
DRAW_BUFFER11 = GL_DRAW_BUFFER11,
DRAW_BUFFER12 = GL_DRAW_BUFFER12,
DRAW_BUFFER13 = GL_DRAW_BUFFER13,
DRAW_BUFFER14 = GL_DRAW_BUFFER14,
DRAW_BUFFER15 = GL_DRAW_BUFFER15,
// returns a single boolean value indicating whether the current edge flag is GL_TRUE or GL_FALSE. The initial value is GL_TRUE. See glEdgeFlag.
EDGE_FLAG = GL_EDGE_FLAG,
// returns a single boolean value indicating whether the edge flag array is enabled. The initial value is GL_FALSE. See glEdgeFlagPointer.
EDGE_FLAG_ARRAY = GL_EDGE_FLAG_ARRAY,
// returns a single value, the name of the buffer object associated with the edge flag array. This buffer object would have been bound to the target GL_ARRAY_BUFFER at the time of the most recent call to glEdgeFlagPointer. If no buffer object was bound to this target, 0 is returned. The initial value is 0. See glBindBuffer.
EDGE_FLAG_ARRAY_BUFFER_BINDING = GL_EDGE_FLAG_ARRAY_BUFFER_BINDING,
// returns one value, the byte offset between consecutive edge flags in the edge flag array. The initial value is 0. See glEdgeFlagPointer.
EDGE_FLAG_ARRAY_STRIDE = GL_EDGE_FLAG_ARRAY_STRIDE,
// returns a single value, the name of the buffer object currently bound to the target GL_ELEMENT_ARRAY_BUFFER. If no buffer object is bound to this target, 0 is returned. The initial value is 0. See glBindBuffer.
ELEMENT_ARRAY_BUFFER_BINDING = GL_ELEMENT_ARRAY_BUFFER_BINDING,
// returns one value, the size of the feedback buffer. See glFeedbackBuffer.
FEEDBACK_BUFFER_SIZE = GL_FEEDBACK_BUFFER_SIZE,
// returns one value, the type of the feedback buffer. See glFeedbackBuffer.
FEEDBACK_BUFFER_TYPE = GL_FEEDBACK_BUFFER_TYPE,
// returns a single boolean value indicating whether fogging is enabled. The initial value is GL_FALSE. See glFog.
FOG = GL_FOG,
// returns a single boolean value indicating whether the fog coordinate array is enabled. The initial value is GL_FALSE. See glFogCoordPointer.
FOG_COORD_ARRAY = GL_FOG_COORD_ARRAY,
// returns a single value, the name of the buffer object associated with the fog coordinate array. This buffer object would have been bound to the target GL_ARRAY_BUFFER at the time of the most recent call to glFogCoordPointer. If no buffer object was bound to this target, 0 is returned. The initial value is 0. See glBindBuffer.
FOG_COORD_ARRAY_BUFFER_BINDING = GL_FOG_COORD_ARRAY_BUFFER_BINDING,
// returns one value, the byte offset between consecutive fog coordinates in the fog coordinate array. The initial value is 0. See glFogCoordPointer.
FOG_COORD_ARRAY_STRIDE = GL_FOG_COORD_ARRAY_STRIDE,
// returns one value, the type of the fog coordinate array. The initial value is GL_FLOAT. See glFogCoordPointer.
FOG_COORD_ARRAY_TYPE = GL_FOG_COORD_ARRAY_TYPE,
// returns one value, a symbolic constant indicating the source of the fog coordinate. The initial value is GL_FRAGMENT_DEPTH. See glFog.
FOG_COORD_SRC = GL_FOG_COORD_SRC,
// returns four values: the red, green, blue, and alpha components of the fog color. Integer values, if requested, are linearly mapped from the internal floating-point representation such that 1.0 returns the most positive representable integer value, and -1.0 returns the most negative representable integer value. The initial value is (0, 0, 0, 0). See glFog.
FOG_COLOR = GL_FOG_COLOR,
// returns one value, the fog density parameter. The initial value is 1. See glFog.
FOG_DENSITY = GL_FOG_DENSITY,
// returns one value, the end factor for the linear fog equation. The initial value is 1. See glFog.
FOG_END = GL_FOG_END,
// returns one value, a symbolic constant indicating the mode of the fog hint. The initial value is GL_DONT_CARE. See glHint.
FOG_HINT = GL_FOG_HINT,
// returns one value, the fog color index. The initial value is 0. See glFog.
FOG_INDEX = GL_FOG_INDEX,
// returns one value, a symbolic constant indicating which fog equation is selected. The initial value is GL_EXP. See glFog.
FOG_MODE = GL_FOG_MODE,
// returns one value, the start factor for the linear fog equation. The initial value is 0. See glFog.
FOG_START = GL_FOG_START,
// returns one value, a symbolic constant indicating the mode of the derivative accuracy hint for fragment shaders. The initial value is GL_DONT_CARE. See glHint.
FRAGMENT_SHADER_DERIVATIVE_HINT = GL_FRAGMENT_SHADER_DERIVATIVE_HINT,
// returns one value, a symbolic constant indicating whether clockwise or counterclockwise polygon winding is treated as front-facing. The initial value is GL_CCW. See glFrontFace.
FRONT_FACE = GL_FRONT_FACE,
// returns one value, a symbolic constant indicating the mode of the mipmap generation filtering hint. The initial value is GL_DONT_CARE. See glHint.
GENERATE_MIPMAP_HINT = GL_GENERATE_MIPMAP_HINT,
// returns one value, the green bias factor used during pixel transfers. The initial value is 0.
GREEN_BIAS = GL_GREEN_BIAS,
// returns one value, the number of green bitplanes in each color buffer.
GREEN_BITS = GL_GREEN_BITS,
// returns one value, the green scale factor used during pixel transfers. The initial value is 1. See glPixelTransfer.
GREEN_SCALE = GL_GREEN_SCALE,
// returns a single boolean value indicating whether histogram is enabled. The initial value is GL_FALSE. See glHistogram.
HISTOGRAM = GL_HISTOGRAM,
// returns a single boolean value indicating whether the color index array is enabled. The initial value is GL_FALSE. See glIndexPointer.
INDEX_ARRAY = GL_INDEX_ARRAY,
// returns a single value, the name of the buffer object associated with the color index array. This buffer object would have been bound to the target GL_ARRAY_BUFFER at the time of the most recent call to glIndexPointer. If no buffer object was bound to this target, 0 is returned. The initial value is 0. See glBindBuffer.
INDEX_ARRAY_BUFFER_BINDING = GL_INDEX_ARRAY_BUFFER_BINDING,
// returns one value, the byte offset between consecutive color indexes in the color index array. The initial value is 0. See glIndexPointer.
INDEX_ARRAY_STRIDE = GL_INDEX_ARRAY_STRIDE,
// returns one value, the data type of indexes in the color index array. The initial value is GL_FLOAT. See glIndexPointer.
INDEX_ARRAY_TYPE = GL_INDEX_ARRAY_TYPE,
// returns one value, the number of bitplanes in each color index buffer.
INDEX_BITS = GL_INDEX_BITS,
// returns one value, the color index used to clear the color index buffers. The initial value is 0. See glClearIndex.
INDEX_CLEAR_VALUE = GL_INDEX_CLEAR_VALUE,
// returns a single boolean value indicating whether a fragment's index values are merged into the framebuffer using a logical operation. The initial value is GL_FALSE. See glLogicOp.
INDEX_LOGIC_OP = GL_INDEX_LOGIC_OP,
// returns a single boolean value indicating whether the GL is in color index mode (GL_TRUE) or RGBA mode (GL_FALSE).
INDEX_MODE = GL_INDEX_MODE,
// returns one value, the offset added to color and stencil indices during pixel transfers. The initial value is 0. See glPixelTransfer.
INDEX_OFFSET = GL_INDEX_OFFSET,
// returns one value, the amount that color and stencil indices are shifted during pixel transfers. The initial value is 0. See glPixelTransfer.
INDEX_SHIFT = GL_INDEX_SHIFT,
// returns one value, a mask indicating which bitplanes of each color index buffer can be written. The initial value is all 1's. See glIndexMask.
INDEX_WRITEMASK = GL_INDEX_WRITEMASK,
// returns a single boolean value indicating whether the specified light is enabled. The initial value is GL_FALSE. See glLight and glLightModel.
LIGHT0 = GL_LIGHT0,
LIGHT1 = GL_LIGHT1,
LIGHT2 = GL_LIGHT2,
LIGHT3 = GL_LIGHT3,
LIGHT4 = GL_LIGHT4,
LIGHT5 = GL_LIGHT5,
LIGHT6 = GL_LIGHT6,
LIGHT7 = GL_LIGHT7,
// returns a single boolean value indicating whether lighting is enabled. The initial value is GL_FALSE. See glLightModel.
LIGHTING = GL_LIGHTING,
// returns four values: the red, green, blue, and alpha components of the ambient intensity of the entire scene. Integer values, if requested, are linearly mapped from the internal floating-point representation such that 1.0 returns the most positive representable integer value, and -1.0 returns the most negative representable integer value. The initial value is (0.2, 0.2, 0.2, 1.0). See glLightModel.
LIGHT_MODEL_AMBIENT = GL_LIGHT_MODEL_AMBIENT,
// returns single enumerated value indicating whether specular reflection calculations are separated from normal lighting computations. The initial value is GL_SINGLE_COLOR.
LIGHT_MODEL_COLOR_CONTROL = GL_LIGHT_MODEL_COLOR_CONTROL,
// returns a single boolean value indicating whether specular reflection calculations treat the viewer as being local to the scene. The initial value is GL_FALSE. See glLightModel.
LIGHT_MODEL_LOCAL_VIEWER = GL_LIGHT_MODEL_LOCAL_VIEWER,
// returns a single boolean value indicating whether separate materials are used to compute lighting for front- and back-facing polygons. The initial value is GL_FALSE. See glLightModel.
LIGHT_MODEL_TWO_SIDE = GL_LIGHT_MODEL_TWO_SIDE,
// returns a single boolean value indicating whether antialiasing of lines is enabled. The initial value is GL_FALSE. See glLineWidth.
LINE_SMOOTH = GL_LINE_SMOOTH,
// returns one value, a symbolic constant indicating the mode of the line antialiasing hint. The initial value is GL_DONT_CARE. See glHint.
LINE_SMOOTH_HINT = GL_LINE_SMOOTH_HINT,
// returns a single boolean value indicating whether stippling of lines is enabled. The initial value is GL_FALSE. See glLineStipple.
LINE_STIPPLE = GL_LINE_STIPPLE,
// returns one value, the 16-bit line stipple pattern. The initial value is all 1's. See glLineStipple.
LINE_STIPPLE_PATTERN = GL_LINE_STIPPLE_PATTERN,
// returns one value, the line stipple repeat factor. The initial value is 1. See glLineStipple.
LINE_STIPPLE_REPEAT = GL_LINE_STIPPLE_REPEAT,
// returns one value, the line width as specified with glLineWidth. The initial value is 1.
LINE_WIDTH = GL_LINE_WIDTH,
// returns one value, the width difference between adjacent supported widths for antialiased lines. See glLineWidth.
LINE_WIDTH_GRANULARITY = GL_LINE_WIDTH_GRANULARITY,
// returns two values: the smallest and largest supported widths for antialiased lines. See glLineWidth.
LINE_WIDTH_RANGE = GL_LINE_WIDTH_RANGE,
// returns one value, the base offset added to all names in arrays presented to glCallLists. The initial value is 0. See glListBase.
LIST_BASE = GL_LIST_BASE,
// returns one value, the name of the display list currently under construction. 0 is returned if no display list is currently under construction. The initial value is 0. See glNewList.
LIST_INDEX = GL_LIST_INDEX,
// returns one value, a symbolic constant indicating the construction mode of the display list currently under construction. The initial value is 0. See glNewList.
LIST_MODE = GL_LIST_MODE,
// returns one value, a symbolic constant indicating the selected logic operation mode. The initial value is GL_COPY. See glLogicOp.
LOGIC_OP_MODE = GL_LOGIC_OP_MODE,
// returns a single boolean value indicating whether 1D evaluation generates colors. The initial value is GL_FALSE. See glMap1.
MAP1_COLOR_4 = GL_MAP1_COLOR_4,
// returns two values: the endpoints of the 1D map's grid domain. The initial value is (0, 1). See glMapGrid.
MAP1_GRID_DOMAIN = GL_MAP1_GRID_DOMAIN,
// returns one value, the number of partitions in the 1D map's grid domain. The initial value is 1. See glMapGrid.
MAP1_GRID_SEGMENTS = GL_MAP1_GRID_SEGMENTS,
// returns a single boolean value indicating whether 1D evaluation generates color indices. The initial value is GL_FALSE. See glMap1.
MAP1_INDEX = GL_MAP1_INDEX,
// returns a single boolean value indicating whether 1D evaluation generates normals. The initial value is GL_FALSE. See glMap1.
MAP1_NORMAL = GL_MAP1_NORMAL,
// returns a single boolean value indicating whether 1D evaluation generates 1D texture coordinates. The initial value is GL_FALSE. See glMap1.
MAP1_TEXTURE_COORD_1 = GL_MAP1_TEXTURE_COORD_1,
// returns a single boolean value indicating whether 1D evaluation generates 2D texture coordinates. The initial value is GL_FALSE. See glMap1.
MAP1_TEXTURE_COORD_2 = GL_MAP1_TEXTURE_COORD_2,
// returns a single boolean value indicating whether 1D evaluation generates 3D texture coordinates. The initial value is GL_FALSE. See glMap1.
MAP1_TEXTURE_COORD_3 = GL_MAP1_TEXTURE_COORD_3,
// returns a single boolean value indicating whether 1D evaluation generates 4D texture coordinates. The initial value is GL_FALSE. See glMap1.
MAP1_TEXTURE_COORD_4 = GL_MAP1_TEXTURE_COORD_4,
// returns a single boolean value indicating whether 1D evaluation generates 3D vertex coordinates. The initial value is GL_FALSE. See glMap1.
MAP1_VERTEX_3 = GL_MAP1_VERTEX_3,
// returns a single boolean value indicating whether 1D evaluation generates 4D vertex coordinates. The initial value is GL_FALSE. See glMap1.
MAP1_VERTEX_4 = GL_MAP1_VERTEX_4,
// returns a single boolean value indicating whether 2D evaluation generates colors. The initial value is GL_FALSE. See glMap2.
MAP2_COLOR_4 = GL_MAP2_COLOR_4,
// returns four values: the endpoints of the 2D map's i and j grid domains. The initial value is (0,1; 0,1). See glMapGrid.
MAP2_GRID_DOMAIN = GL_MAP2_GRID_DOMAIN,
// returns two values: the number of partitions in the 2D map's i and j grid domains. The initial value is (1,1). See glMapGrid.
MAP2_GRID_SEGMENTS = GL_MAP2_GRID_SEGMENTS,
// returns a single boolean value indicating whether 2D evaluation generates color indices. The initial value is GL_FALSE. See glMap2.
MAP2_INDEX = GL_MAP2_INDEX,
// returns a single boolean value indicating whether 2D evaluation generates normals. The initial value is GL_FALSE. See glMap2.
MAP2_NORMAL = GL_MAP2_NORMAL,
// returns a single boolean value indicating whether 2D evaluation generates 1D texture coordinates. The initial value is GL_FALSE. See glMap2.
MAP2_TEXTURE_COORD_1 = GL_MAP2_TEXTURE_COORD_1,
// returns a single boolean value indicating whether 2D evaluation generates 2D texture coordinates. The initial value is GL_FALSE. See glMap2.
MAP2_TEXTURE_COORD_2 = GL_MAP2_TEXTURE_COORD_2,
// returns a single boolean value indicating whether 2D evaluation generates 3D texture coordinates. The initial value is GL_FALSE. See glMap2.
MAP2_TEXTURE_COORD_3 = GL_MAP2_TEXTURE_COORD_3,
// returns a single boolean value indicating whether 2D evaluation generates 4D texture coordinates. The initial value is GL_FALSE. See glMap2.
MAP2_TEXTURE_COORD_4 = GL_MAP2_TEXTURE_COORD_4,
// returns a single boolean value indicating whether 2D evaluation generates 3D vertex coordinates. The initial value is GL_FALSE. See glMap2.
MAP2_VERTEX_3 = GL_MAP2_VERTEX_3,
// returns a single boolean value indicating whether 2D evaluation generates 4D vertex coordinates. The initial value is GL_FALSE. See glMap2.
MAP2_VERTEX_4 = GL_MAP2_VERTEX_4,
// returns a single boolean value indicating if colors and color indices are to be replaced by table lookup during pixel transfers. The initial value is GL_FALSE. See glPixelTransfer.
MAP_COLOR = GL_MAP_COLOR,
// returns a single boolean value indicating if stencil indices are to be replaced by table lookup during pixel transfers. The initial value is GL_FALSE. See glPixelTransfer.
MAP_STENCIL = GL_MAP_STENCIL,
// returns one value, a symbolic constant indicating which matrix stack is currently the target of all matrix operations. The initial value is GL_MODELVIEW. See glMatrixMode.
MATRIX_MODE = GL_MATRIX_MODE,
// returns one value, a rough estimate of the largest 3D texture that the GL can handle. The value must be at least 16. If the GL version is 1.2 or greater, use GL_PROXY_TEXTURE_3D to determine if a texture is too large. See glTexImage3D.
MAX_3D_TEXTURE_SIZE = GL_MAX_3D_TEXTURE_SIZE,
// returns one value indicating the maximum supported depth of the client attribute stack. See glPushClientAttrib.
MAX_CLIENT_ATTRIB_STACK_DEPTH = GL_MAX_CLIENT_ATTRIB_STACK_DEPTH,
// returns one value, the maximum supported depth of the attribute stack. The value must be at least 16. See glPushAttrib.
MAX_ATTRIB_STACK_DEPTH = GL_MAX_ATTRIB_STACK_DEPTH,
// returns one value, the maximum number of application-defined clipping planes. The value must be at least 6. See glClipPlane.
MAX_CLIP_PLANES = GL_MAX_CLIP_PLANES,
// returns one value, the maximum supported depth of the color matrix stack. The value must be at least 2. See glPushMatrix.
MAX_COLOR_MATRIX_STACK_DEPTH = GL_MAX_COLOR_MATRIX_STACK_DEPTH,
// returns one value, the maximum supported texture image units that can be used to access texture maps from the vertex shader and the fragment processor combined. If both the vertex shader and the fragment processing stage access the same texture image unit, then that counts as using two texture image units against this limit. The value must be at least 2. See glActiveTexture.
MAX_COMBINED_TEXTURE_IMAGE_UNITS = GL_MAX_COMBINED_TEXTURE_IMAGE_UNITS,
// returns one value. The value gives a rough estimate of the largest cube-map texture that the GL can handle. The value must be at least 16. If the GL version is 1.3 or greater, use GL_PROXY_TEXTURE_CUBE_MAP to determine if a texture is too large. See glTexImage2D.
MAX_CUBE_MAP_TEXTURE_SIZE = GL_MAX_CUBE_MAP_TEXTURE_SIZE,
// returns one value, the maximum number of simultaneous output colors allowed from a fragment shader using the gl_FragData built-in array. The value must be at least 1. See glDrawBuffers.
MAX_DRAW_BUFFERS = GL_MAX_DRAW_BUFFERS,
// returns one value, the recommended maximum number of vertex array indices. See glDrawRangeElements.
MAX_ELEMENTS_INDICES = GL_MAX_ELEMENTS_INDICES,
// returns one value, the recommended maximum number of vertex array vertices. See glDrawRangeElements.
MAX_ELEMENTS_VERTICES = GL_MAX_ELEMENTS_VERTICES,
// returns one value, the maximum equation order supported by 1D and 2D evaluators. The value must be at least 8. See glMap1 and glMap2.
MAX_EVAL_ORDER = GL_MAX_EVAL_ORDER,
// returns one value, the maximum number of individual floating-point, integer, or boolean values that can be held in uniform variable storage for a fragment shader. The value must be at least 64. See glUniform.
MAX_FRAGMENT_UNIFORM_COMPONENTS = GL_MAX_FRAGMENT_UNIFORM_COMPONENTS,
// returns one value, the maximum number of lights. The value must be at least 8. See glLight.
MAX_LIGHTS = GL_MAX_LIGHTS,
// returns one value, the maximum recursion depth allowed during display-list traversal. The value must be at least 64. See glCallList.
MAX_LIST_NESTING = GL_MAX_LIST_NESTING,
// returns one value, the maximum supported depth of the modelview matrix stack. The value must be at least 32. See glPushMatrix.
MAX_MODELVIEW_STACK_DEPTH = GL_MAX_MODELVIEW_STACK_DEPTH,
// returns one value, the maximum supported depth of the selection name stack. The value must be at least 64. See glPushName.
MAX_NAME_STACK_DEPTH = GL_MAX_NAME_STACK_DEPTH,
// returns one value, the maximum supported size of a glPixelMap lookup table. The value must be at least 32. See glPixelMap.
MAX_PIXEL_MAP_TABLE = GL_MAX_PIXEL_MAP_TABLE,
// returns one value, the maximum supported depth of the projection matrix stack. The value must be at least 2. See glPushMatrix.
MAX_PROJECTION_STACK_DEPTH = GL_MAX_PROJECTION_STACK_DEPTH,
// returns one value, the maximum number of texture coordinate sets available to vertex and fragment shaders. The value must be at least 2. See glActiveTexture and glClientActiveTexture.
MAX_TEXTURE_COORDS = GL_MAX_TEXTURE_COORDS,
// returns one value, the maximum supported texture image units that can be used to access texture maps from the fragment shader. The value must be at least 2. See glActiveTexture.
MAX_TEXTURE_IMAGE_UNITS = GL_MAX_TEXTURE_IMAGE_UNITS,
// returns one value, the maximum, absolute value of the texture level-of-detail bias. The value must be at least 4.
MAX_TEXTURE_LOD_BIAS = GL_MAX_TEXTURE_LOD_BIAS,
// returns one value. The value gives a rough estimate of the largest texture that the GL can handle. The value must be at least 64. If the GL version is 1.1 or greater, use GL_PROXY_TEXTURE_1D or GL_PROXY_TEXTURE_2D to determine if a texture is too large. See glTexImage1D and glTexImage2D.
MAX_TEXTURE_SIZE = GL_MAX_TEXTURE_SIZE,
// returns one value, the maximum supported depth of the texture matrix stack. The value must be at least 2. See glPushMatrix.
MAX_TEXTURE_STACK_DEPTH = GL_MAX_TEXTURE_STACK_DEPTH,
// returns a single value indicating the number of conventional texture units supported. Each conventional texture unit includes both a texture coordinate set and a texture image unit. Conventional texture units may be used for fixed-function (non-shader) rendering. The value must be at least 2. Additional texture coordinate sets and texture image units may be accessed from vertex and fragment shaders. See glActiveTexture and glClientActiveTexture.
MAX_TEXTURE_UNITS = GL_MAX_TEXTURE_UNITS,
// returns one value, the maximum number of interpolators available for processing varying variables used by vertex and fragment shaders. This value represents the number of individual floating-point values that can be interpolated; varying variables declared as vectors, matrices, and arrays will all consume multiple interpolators. The value must be at least 32.
MAX_VARYING_FLOATS = GL_MAX_VARYING_FLOATS,
// returns one value, the maximum number of 4-component generic vertex attributes accessible to a vertex shader. The value must be at least 16. See glVertexAttrib.
MAX_VERTEX_ATTRIBS = GL_MAX_VERTEX_ATTRIBS,
// returns one value, the maximum supported texture image units that can be used to access texture maps from the vertex shader. The value may be 0. See glActiveTexture.
MAX_VERTEX_TEXTURE_IMAGE_UNITS = GL_MAX_VERTEX_TEXTURE_IMAGE_UNITS,
// returns one value, the maximum number of individual floating-point, integer, or boolean values that can be held in uniform variable storage for a vertex shader. The value must be at least 512. See glUniform.
MAX_VERTEX_UNIFORM_COMPONENTS = GL_MAX_VERTEX_UNIFORM_COMPONENTS,
// returns two values: the maximum supported width and height of the viewport. These must be at least as large as the visible dimensions of the display being rendered to. See glViewport.
MAX_VIEWPORT_DIMS = GL_MAX_VIEWPORT_DIMS,
// returns a single boolean value indicating whether pixel minmax values are computed. The initial value is GL_FALSE. See glMinmax.
MINMAX = GL_MINMAX,
// returns sixteen values: the modelview matrix on the top of the modelview matrix stack. Initially this matrix is the identity matrix. See glPushMatrix.
MODELVIEW_MATRIX = GL_MODELVIEW_MATRIX,
// returns one value, the number of matrices on the modelview matrix stack. The initial value is 1. See glPushMatrix.
MODELVIEW_STACK_DEPTH = GL_MODELVIEW_STACK_DEPTH,
// returns one value, the number of names on the selection name stack. The initial value is 0. See glPushName.
NAME_STACK_DEPTH = GL_NAME_STACK_DEPTH,
// returns a single boolean value, indicating whether the normal array is enabled. The initial value is GL_FALSE. See glNormalPointer.
NORMAL_ARRAY = GL_NORMAL_ARRAY,
// returns a single value, the name of the buffer object associated with the normal array. This buffer object would have been bound to the target GL_ARRAY_BUFFER at the time of the most recent call to glNormalPointer. If no buffer object was bound to this target, 0 is returned. The initial value is 0. See glBindBuffer.
NORMAL_ARRAY_BUFFER_BINDING = GL_NORMAL_ARRAY_BUFFER_BINDING,
// returns one value, the byte offset between consecutive normals in the normal array. The initial value is 0. See glNormalPointer.
NORMAL_ARRAY_STRIDE = GL_NORMAL_ARRAY_STRIDE,
// returns one value, the data type of each coordinate in the normal array. The initial value is GL_FLOAT. See glNormalPointer.
NORMAL_ARRAY_TYPE = GL_NORMAL_ARRAY_TYPE,
// returns a single boolean value indicating whether normals are automatically scaled to unit length after they have been transformed to eye coordinates. The initial value is GL_FALSE. See glNormal.
NORMALIZE = GL_NORMALIZE,
// returns a single integer value indicating the number of available compressed texture formats. The minimum value is 0. See glCompressedTexImage2D.
NUM_COMPRESSED_TEXTURE_FORMATS = GL_NUM_COMPRESSED_TEXTURE_FORMATS,
// returns one value, the byte alignment used for writing pixel data to memory. The initial value is 4. See glPixelStore.
PACK_ALIGNMENT = GL_PACK_ALIGNMENT,
// returns one value, the image height used for writing pixel data to memory. The initial value is 0. See glPixelStore.
PACK_IMAGE_HEIGHT = GL_PACK_IMAGE_HEIGHT,
// returns a single boolean value indicating whether single-bit pixels being written to memory are written first to the least significant bit of each unsigned byte. The initial value is GL_FALSE. See glPixelStore.
PACK_LSB_FIRST = GL_PACK_LSB_FIRST,
// returns one value, the row length used for writing pixel data to memory. The initial value is 0. See glPixelStore.
PACK_ROW_LENGTH = GL_PACK_ROW_LENGTH,
// returns one value, the number of pixel images skipped before the first pixel is written into memory. The initial value is 0. See glPixelStore.
PACK_SKIP_IMAGES = GL_PACK_SKIP_IMAGES,
// returns one value, the number of pixel locations skipped before the first pixel is written into memory. The initial value is 0. See glPixelStore.
PACK_SKIP_PIXELS = GL_PACK_SKIP_PIXELS,
// returns one value, the number of rows of pixel locations skipped before the first pixel is written into memory. The initial value is 0. See glPixelStore.
PACK_SKIP_ROWS = GL_PACK_SKIP_ROWS,
// returns a single boolean value indicating whether the bytes of two-byte and four-byte pixel indices and components are swapped before being written to memory. The initial value is GL_FALSE. See glPixelStore.
PACK_SWAP_BYTES = GL_PACK_SWAP_BYTES,
// returns one value, a symbolic constant indicating the mode of the perspective correction hint. The initial value is GL_DONT_CARE. See glHint.
PERSPECTIVE_CORRECTION_HINT = GL_PERSPECTIVE_CORRECTION_HINT,
// returns one value, the size of the alpha-to-alpha pixel translation table. The initial value is 1. See glPixelMap.
PIXEL_MAP_A_TO_A_SIZE = GL_PIXEL_MAP_A_TO_A_SIZE,
// returns one value, the size of the blue-to-blue pixel translation table. The initial value is 1. See glPixelMap.
PIXEL_MAP_B_TO_B_SIZE = GL_PIXEL_MAP_B_TO_B_SIZE,
// returns one value, the size of the green-to-green pixel translation table. The initial value is 1. See glPixelMap.
PIXEL_MAP_G_TO_G_SIZE = GL_PIXEL_MAP_G_TO_G_SIZE,
// returns one value, the size of the index-to-alpha pixel translation table. The initial value is 1. See glPixelMap.
PIXEL_MAP_I_TO_A_SIZE = GL_PIXEL_MAP_I_TO_A_SIZE,
// returns one value, the size of the index-to-blue pixel translation table. The initial value is 1. See glPixelMap.
PIXEL_MAP_I_TO_B_SIZE = GL_PIXEL_MAP_I_TO_B_SIZE,
// returns one value, the size of the index-to-green pixel translation table. The initial value is 1. See glPixelMap.
PIXEL_MAP_I_TO_G_SIZE = GL_PIXEL_MAP_I_TO_G_SIZE,
// returns one value, the size of the index-to-index pixel translation table. The initial value is 1. See glPixelMap.
PIXEL_MAP_I_TO_I_SIZE = GL_PIXEL_MAP_I_TO_I_SIZE,
// returns one value, the size of the index-to-red pixel translation table. The initial value is 1. See glPixelMap.
PIXEL_MAP_I_TO_R_SIZE = GL_PIXEL_MAP_I_TO_R_SIZE,
// returns one value, the size of the red-to-red pixel translation table. The initial value is 1. See glPixelMap.
PIXEL_MAP_R_TO_R_SIZE = GL_PIXEL_MAP_R_TO_R_SIZE,
// returns one value, the size of the stencil-to-stencil pixel translation table. The initial value is 1. See glPixelMap.
PIXEL_MAP_S_TO_S_SIZE = GL_PIXEL_MAP_S_TO_S_SIZE,
// returns a single value, the name of the buffer object currently bound to the target GL_PIXEL_PACK_BUFFER. If no buffer object is bound to this target, 0 is returned. The initial value is 0. See glBindBuffer.
PIXEL_PACK_BUFFER_BINDING = GL_PIXEL_PACK_BUFFER_BINDING,
// returns a single value, the name of the buffer object currently bound to the target GL_PIXEL_UNPACK_BUFFER. If no buffer object is bound to this target, 0 is returned. The initial value is 0. See glBindBuffer.
PIXEL_UNPACK_BUFFER_BINDING = GL_PIXEL_UNPACK_BUFFER_BINDING,
// returns three values, the coefficients for computing the attenuation value for points. See glPointParameter.
POINT_DISTANCE_ATTENUATION = GL_POINT_DISTANCE_ATTENUATION,
// returns one value, the point size threshold for determining the point size. See glPointParameter.
POINT_FADE_THRESHOLD_SIZE = GL_POINT_FADE_THRESHOLD_SIZE,
// returns one value, the point size as specified by glPointSize. The initial value is 1.
POINT_SIZE = GL_POINT_SIZE,
// returns one value, the size difference between adjacent supported sizes for antialiased points. See glPointSize.
POINT_SIZE_GRANULARITY = GL_POINT_SIZE_GRANULARITY,
// returns one value, the upper bound for the attenuated point sizes. The initial value is 0.0. See glPointParameter.
POINT_SIZE_MAX = GL_POINT_SIZE_MAX,
// returns one value, the lower bound for the attenuated point sizes. The initial value is 1.0. See glPointParameter.
POINT_SIZE_MIN = GL_POINT_SIZE_MIN,
// returns two values: the smallest and largest supported sizes for antialiased points. The smallest size must be at most 1, and the largest size must be at least 1. See glPointSize.
POINT_SIZE_RANGE = GL_POINT_SIZE_RANGE,
// returns a single boolean value indicating whether antialiasing of points is enabled. The initial value is GL_FALSE. See glPointSize.
POINT_SMOOTH = GL_POINT_SMOOTH,
// returns one value, a symbolic constant indicating the mode of the point antialiasing hint. The initial value is GL_DONT_CARE. See glHint.
POINT_SMOOTH_HINT = GL_POINT_SMOOTH_HINT,
// returns a single boolean value indicating whether point sprite is enabled. The initial value is GL_FALSE.
POINT_SPRITE = GL_POINT_SPRITE,
// returns two values: symbolic constants indicating whether front-facing and back-facing polygons are rasterized as points, lines, or filled polygons. The initial value is GL_FILL. See glPolygonMode.
POLYGON_MODE = GL_POLYGON_MODE,
// returns one value, the scaling factor used to determine the variable offset that is added to the depth value of each fragment generated when a polygon is rasterized. The initial value is 0. See glPolygonOffset.
POLYGON_OFFSET_FACTOR = GL_POLYGON_OFFSET_FACTOR,
// returns one value. This value is multiplied by an implementation-specific value and then added to the depth value of each fragment generated when a polygon is rasterized. The initial value is 0. See glPolygonOffset.
POLYGON_OFFSET_UNITS = GL_POLYGON_OFFSET_UNITS,
// returns a single boolean value indicating whether polygon offset is enabled for polygons in fill mode. The initial value is GL_FALSE. See glPolygonOffset.
POLYGON_OFFSET_FILL = GL_POLYGON_OFFSET_FILL,
// returns a single boolean value indicating whether polygon offset is enabled for polygons in line mode. The initial value is GL_FALSE. See glPolygonOffset.
POLYGON_OFFSET_LINE = GL_POLYGON_OFFSET_LINE,
// returns a single boolean value indicating whether polygon offset is enabled for polygons in point mode. The initial value is GL_FALSE. See glPolygonOffset.
POLYGON_OFFSET_POINT = GL_POLYGON_OFFSET_POINT,
// returns a single boolean value indicating whether antialiasing of polygons is enabled. The initial value is GL_FALSE. See glPolygonMode.
POLYGON_SMOOTH = GL_POLYGON_SMOOTH,
// returns one value, a symbolic constant indicating the mode of the polygon antialiasing hint. The initial value is GL_DONT_CARE. See glHint.
POLYGON_SMOOTH_HINT = GL_POLYGON_SMOOTH_HINT,
// returns a single boolean value indicating whether polygon stippling is enabled. The initial value is GL_FALSE. See glPolygonStipple.
POLYGON_STIPPLE = GL_POLYGON_STIPPLE,
// returns a single boolean value indicating whether post color matrix transformation lookup is enabled. The initial value is GL_FALSE. See glColorTable.
POST_COLOR_MATRIX_COLOR_TABLE = GL_POST_COLOR_MATRIX_COLOR_TABLE,
// returns one value, the red bias factor applied to RGBA fragments after color matrix transformations. The initial value is 0. See glPixelTransfer.
POST_COLOR_MATRIX_RED_BIAS = GL_POST_COLOR_MATRIX_RED_BIAS,
// returns one value, the green bias factor applied to RGBA fragments after color matrix transformations. The initial value is 0. See glPixelTransfer
POST_COLOR_MATRIX_GREEN_BIAS = GL_POST_COLOR_MATRIX_GREEN_BIAS,
// returns one value, the blue bias factor applied to RGBA fragments after color matrix transformations. The initial value is 0. See glPixelTransfer.
POST_COLOR_MATRIX_BLUE_BIAS = GL_POST_COLOR_MATRIX_BLUE_BIAS,
// returns one value, the alpha bias factor applied to RGBA fragments after color matrix transformations. The initial value is 0. See glPixelTransfer.
POST_COLOR_MATRIX_ALPHA_BIAS = GL_POST_COLOR_MATRIX_ALPHA_BIAS,
// returns one value, the red scale factor applied to RGBA fragments after color matrix transformations. The initial value is 1. See glPixelTransfer.
POST_COLOR_MATRIX_RED_SCALE = GL_POST_COLOR_MATRIX_RED_SCALE,
// returns one value, the green scale factor applied to RGBA fragments after color matrix transformations. The initial value is 1. See glPixelTransfer.
POST_COLOR_MATRIX_GREEN_SCALE = GL_POST_COLOR_MATRIX_GREEN_SCALE,
// returns one value, the blue scale factor applied to RGBA fragments after color matrix transformations. The initial value is 1. See glPixelTransfer.
POST_COLOR_MATRIX_BLUE_SCALE = GL_POST_COLOR_MATRIX_BLUE_SCALE,
// returns one value, the alpha scale factor applied to RGBA fragments after color matrix transformations. The initial value is 1. See glPixelTransfer.
POST_COLOR_MATRIX_ALPHA_SCALE = GL_POST_COLOR_MATRIX_ALPHA_SCALE,
// returns a single boolean value indicating whether post convolution lookup is enabled. The initial value is GL_FALSE. See glColorTable.
POST_CONVOLUTION_COLOR_TABLE = GL_POST_CONVOLUTION_COLOR_TABLE,
// returns one value, the red bias factor applied to RGBA fragments after convolution. The initial value is 0. See glPixelTransfer.
POST_CONVOLUTION_RED_BIAS = GL_POST_CONVOLUTION_RED_BIAS,
// returns one value, the green bias factor applied to RGBA fragments after convolution. The initial value is 0. See glPixelTransfer.
POST_CONVOLUTION_GREEN_BIAS = GL_POST_CONVOLUTION_GREEN_BIAS,
// returns one value, the blue bias factor applied to RGBA fragments after convolution. The initial value is 0. See glPixelTransfer.
POST_CONVOLUTION_BLUE_BIAS = GL_POST_CONVOLUTION_BLUE_BIAS,
// returns one value, the alpha bias factor applied to RGBA fragments after convolution. The initial value is 0. See glPixelTransfer.
POST_CONVOLUTION_ALPHA_BIAS = GL_POST_CONVOLUTION_ALPHA_BIAS,
// returns one value, the red scale factor applied to RGBA fragments after convolution. The initial value is 1. See glPixelTransfer.
POST_CONVOLUTION_RED_SCALE = GL_POST_CONVOLUTION_RED_SCALE,
// returns one value, the green scale factor applied to RGBA fragments after convolution. The initial value is 1. See glPixelTransfer.
POST_CONVOLUTION_GREEN_SCALE = GL_POST_CONVOLUTION_GREEN_SCALE,
// returns one value, the blue scale factor applied to RGBA fragments after convolution. The initial value is 1. See glPixelTransfer.
POST_CONVOLUTION_BLUE_SCALE = GL_POST_CONVOLUTION_BLUE_SCALE,
// returns one value, the alpha scale factor applied to RGBA fragments after convolution. The initial value is 1. See glPixelTransfer.
POST_CONVOLUTION_ALPHA_SCALE = GL_POST_CONVOLUTION_ALPHA_SCALE,
// returns sixteen values: the projection matrix on the top of the projection matrix stack. Initially this matrix is the identity matrix. See glPushMatrix.
PROJECTION_MATRIX = GL_PROJECTION_MATRIX,
// returns one value, the number of matrices on the projection matrix stack. The initial value is 1. See glPushMatrix.
PROJECTION_STACK_DEPTH = GL_PROJECTION_STACK_DEPTH,
// returns one value, a symbolic constant indicating which color buffer is selected for reading. The initial value is GL_BACK if there is a back buffer, otherwise it is GL_FRONT. See glReadPixels and glAccum.
READ_BUFFER = GL_READ_BUFFER,
// returns one value, the red bias factor used during pixel transfers. The initial value is 0.
RED_BIAS = GL_RED_BIAS,
// returns one value, the number of red bitplanes in each color buffer.
RED_BITS = GL_RED_BITS,
// returns one value, the red scale factor used during pixel transfers. The initial value is 1. See glPixelTransfer.
RED_SCALE = GL_RED_SCALE,
// returns one value, a symbolic constant indicating whether the GL is in render, select, or feedback mode. The initial value is GL_RENDER. See glRenderMode.
RENDER_MODE = GL_RENDER_MODE,
// returns single boolean value indicating whether normal rescaling is enabled. See glEnable.
RESCALE_NORMAL = GL_RESCALE_NORMAL,
// returns a single boolean value indicating whether the GL is in RGBA mode (true) or color index mode (false). See glColor.
RGBA_MODE = GL_RGBA_MODE,
// returns a single integer value indicating the number of sample buffers associated with the framebuffer. See glSampleCoverage.
SAMPLE_BUFFERS = GL_SAMPLE_BUFFERS,
// returns a single positive floating-point value indicating the current sample coverage value. See glSampleCoverage.
SAMPLE_COVERAGE_VALUE = GL_SAMPLE_COVERAGE_VALUE,
// returns a single boolean value indicating if the temporary coverage value should be inverted. See glSampleCoverage.
SAMPLE_COVERAGE_INVERT = GL_SAMPLE_COVERAGE_INVERT,
// returns a single integer value indicating the coverage mask size. See glSampleCoverage.
SAMPLES = GL_SAMPLES,
// returns four values: the x and y window coordinates of the scissor box, followed by its width and height. Initially the x and y window coordinates are both 0 and the width and height are set to the size of the window. See glScissor.
SCISSOR_BOX = GL_SCISSOR_BOX,
// returns a single boolean value indicating whether scissoring is enabled. The initial value is GL_FALSE. See glScissor.
SCISSOR_TEST = GL_SCISSOR_TEST,
// returns a single boolean value indicating whether the secondary color array is enabled. The initial value is GL_FALSE. See glSecondaryColorPointer.
SECONDARY_COLOR_ARRAY = GL_SECONDARY_COLOR_ARRAY,
// returns a single value, the name of the buffer object associated with the secondary color array. This buffer object would have been bound to the target GL_ARRAY_BUFFER at the time of the most recent call to glSecondaryColorPointer. If no buffer object was bound to this target, 0 is returned. The initial value is 0. See glBindBuffer.
SECONDARY_COLOR_ARRAY_BUFFER_BINDING = GL_SECONDARY_COLOR_ARRAY_BUFFER_BINDING,
// returns one value, the number of components per color in the secondary color array. The initial value is 3. See glSecondaryColorPointer.
SECONDARY_COLOR_ARRAY_SIZE = GL_SECONDARY_COLOR_ARRAY_SIZE,
// returns one value, the byte offset between consecutive colors in the secondary color array. The initial value is 0. See glSecondaryColorPointer.
SECONDARY_COLOR_ARRAY_STRIDE = GL_SECONDARY_COLOR_ARRAY_STRIDE,
// returns one value, the data type of each component in the secondary color array. The initial value is GL_FLOAT. See glSecondaryColorPointer.
SECONDARY_COLOR_ARRAY_TYPE = GL_SECONDARY_COLOR_ARRAY_TYPE,
// return one value, the size of the selection buffer. See glSelectBuffer.
SELECTION_BUFFER_SIZE = GL_SELECTION_BUFFER_SIZE,
// returns a single boolean value indicating whether 2D separable convolution is enabled. The initial value is GL_FALSE. See glSeparableFilter2D.
SEPARABLE_2D = GL_SEPARABLE_2D,
// returns one value, a symbolic constant indicating whether the shading mode is flat or smooth. The initial value is GL_SMOOTH. See glShadeModel.
SHADE_MODEL = GL_SHADE_MODEL,
// returns two values, the smallest and largest supported widths for antialiased lines. See glLineWidth.
SMOOTH_LINE_WIDTH_RANGE = GL_SMOOTH_LINE_WIDTH_RANGE,
// returns one value, the granularity of widths for antialiased lines. See glLineWidth.
SMOOTH_LINE_WIDTH_GRANULARITY = GL_SMOOTH_LINE_WIDTH_GRANULARITY,
// returns two values, the smallest and largest supported widths for antialiased points. See glPointSize.
SMOOTH_POINT_SIZE_RANGE = GL_SMOOTH_POINT_SIZE_RANGE,
// returns one value, the granularity of sizes for antialiased points. See glPointSize.
SMOOTH_POINT_SIZE_GRANULARITY = GL_SMOOTH_POINT_SIZE_GRANULARITY,
// returns one value, a symbolic constant indicating what action is taken for back-facing polygons when the stencil test fails. The initial value is GL_KEEP. See glStencilOpSeparate.
STENCIL_BACK_FAIL = GL_STENCIL_BACK_FAIL,
// returns one value, a symbolic constant indicating what function is used for back-facing polygons to compare the stencil reference value with the stencil buffer value. The initial value is GL_ALWAYS. See glStencilFuncSeparate.
STENCIL_BACK_FUNC = GL_STENCIL_BACK_FUNC,
// returns one value, a symbolic constant indicating what action is taken for back-facing polygons when the stencil test passes, but the depth test fails. The initial value is GL_KEEP. See glStencilOpSeparate.
STENCIL_BACK_PASS_DEPTH_FAIL = GL_STENCIL_BACK_PASS_DEPTH_FAIL,
// returns one value, a symbolic constant indicating what action is taken for back-facing polygons when the stencil test passes and the depth test passes. The initial value is GL_KEEP. See glStencilOpSeparate.
STENCIL_BACK_PASS_DEPTH_PASS = GL_STENCIL_BACK_PASS_DEPTH_PASS,
// returns one value, the reference value that is compared with the contents of the stencil buffer for back-facing polygons. The initial value is 0. See glStencilFuncSeparate.
STENCIL_BACK_REF = GL_STENCIL_BACK_REF,
// returns one value, the mask that is used for back-facing polygons to mask both the stencil reference value and the stencil buffer value before they are compared. The initial value is all 1's. See glStencilFuncSeparate.
STENCIL_BACK_VALUE_MASK = GL_STENCIL_BACK_VALUE_MASK,
// returns one value, the mask that controls writing of the stencil bitplanes for back-facing polygons. The initial value is all 1's. See glStencilMaskSeparate.
STENCIL_BACK_WRITEMASK = GL_STENCIL_BACK_WRITEMASK,
// returns one value, the number of bitplanes in the stencil buffer.
STENCIL_BITS = GL_STENCIL_BITS,
// returns one value, the index to which the stencil bitplanes are cleared. The initial value is 0. See glClearStencil.
STENCIL_CLEAR_VALUE = GL_STENCIL_CLEAR_VALUE,
// returns one value, a symbolic constant indicating what action is taken when the stencil test fails. The initial value is GL_KEEP. See glStencilOp. If the GL version is 2.0 or greater, this stencil state only affects non-polygons and front-facing polygons. Back-facing polygons use separate stencil state. See glStencilOpSeparate.
STENCIL_FAIL = GL_STENCIL_FAIL,
// returns one value, a symbolic constant indicating what function is used to compare the stencil reference value with the stencil buffer value. The initial value is GL_ALWAYS. See glStencilFunc. If the GL version is 2.0 or greater, this stencil state only affects non-polygons and front-facing polygons. Back-facing polygons use separate stencil state. See glStencilFuncSeparate.
STENCIL_FUNC = GL_STENCIL_FUNC,
// returns one value, a symbolic constant indicating what action is taken when the stencil test passes, but the depth test fails. The initial value is GL_KEEP. See glStencilOp. If the GL version is 2.0 or greater, this stencil state only affects non-polygons and front-facing polygons. Back-facing polygons use separate stencil state. See glStencilOpSeparate.
STENCIL_PASS_DEPTH_FAIL = GL_STENCIL_PASS_DEPTH_FAIL,
// returns one value, a symbolic constant indicating what action is taken when the stencil test passes and the depth test passes. The initial value is GL_KEEP. See glStencilOp. If the GL version is 2.0 or greater, this stencil state only affects non-polygons and front-facing polygons. Back-facing polygons use separate stencil state. See glStencilOpSeparate.
STENCIL_PASS_DEPTH_PASS = GL_STENCIL_PASS_DEPTH_PASS,
// returns one value, the reference value that is compared with the contents of the stencil buffer. The initial value is 0. See glStencilFunc. If the GL version is 2.0 or greater, this stencil state only affects non-polygons and front-facing polygons. Back-facing polygons use separate stencil state. See glStencilFuncSeparate.
STENCIL_REF = GL_STENCIL_REF,
// returns a single boolean value indicating whether stencil testing of fragments is enabled. The initial value is GL_FALSE. See glStencilFunc and glStencilOp.
STENCIL_TEST = GL_STENCIL_TEST,
// returns one value, the mask that is used to mask both the stencil reference value and the stencil buffer value before they are compared. The initial value is all 1's. See glStencilFunc. If the GL version is 2.0 or greater, this stencil state only affects non-polygons and front-facing polygons. Back-facing polygons use separate stencil state. See glStencilFuncSeparate.
STENCIL_VALUE_MASK = GL_STENCIL_VALUE_MASK,
// returns one value, the mask that controls writing of the stencil bitplanes. The initial value is all 1's. See glStencilMask. If the GL version is 2.0 or greater, this stencil state only affects non-polygons and front-facing polygons. Back-facing polygons use separate stencil state. See glStencilMaskSeparate.
STENCIL_WRITEMASK = GL_STENCIL_WRITEMASK,
// returns a single boolean value indicating whether stereo buffers (left and right) are supported.
STEREO = GL_STEREO,
// returns one value, an estimate of the number of bits of subpixel resolution that are used to position rasterized geometry in window coordinates. The value must be at least 4.
SUBPIXEL_BITS = GL_SUBPIXEL_BITS,
// returns a single boolean value indicating whether 1D texture mapping is enabled. The initial value is GL_FALSE. See glTexImage1D.
TEXTURE_1D = GL_TEXTURE_1D,
// returns a single value, the name of the texture currently bound to the target GL_TEXTURE_1D. The initial value is 0. See glBindTexture.
TEXTURE_BINDING_1D = GL_TEXTURE_BINDING_1D,
// returns a single boolean value indicating whether 2D texture mapping is enabled. The initial value is GL_FALSE. See glTexImage2D.
TEXTURE_2D = GL_TEXTURE_2D,
// returns a single value, the name of the texture currently bound to the target GL_TEXTURE_2D. The initial value is 0. See glBindTexture.
TEXTURE_BINDING_2D = GL_TEXTURE_BINDING_2D,
// returns a single boolean value indicating whether 3D texture mapping is enabled. The initial value is GL_FALSE. See glTexImage3D.
TEXTURE_3D = GL_TEXTURE_3D,
// returns a single value, the name of the texture currently bound to the target GL_TEXTURE_3D. The initial value is 0. See glBindTexture.
TEXTURE_BINDING_3D = GL_TEXTURE_BINDING_3D,
// returns a single value, the name of the texture currently bound to the target GL_TEXTURE_CUBE_MAP. The initial value is 0. See glBindTexture.
TEXTURE_BINDING_CUBE_MAP = GL_TEXTURE_BINDING_CUBE_MAP,
// returns a single value indicating the mode of the texture compression hint. The initial value is GL_DONT_CARE.
TEXTURE_COMPRESSION_HINT = GL_TEXTURE_COMPRESSION_HINT,
// returns a single boolean value indicating whether the texture coordinate array is enabled. The initial value is GL_FALSE. See glTexCoordPointer.
TEXTURE_COORD_ARRAY = GL_TEXTURE_COORD_ARRAY,
// returns a single value, the name of the buffer object associated with the texture coordinate array. This buffer object would have been bound to the target GL_ARRAY_BUFFER at the time of the most recent call to glTexCoordPointer. If no buffer object was bound to this target, 0 is returned. The initial value is 0. See glBindBuffer.
TEXTURE_COORD_ARRAY_BUFFER_BINDING = GL_TEXTURE_COORD_ARRAY_BUFFER_BINDING,
// returns one value, the number of coordinates per element in the texture coordinate array. The initial value is 4. See glTexCoordPointer.
TEXTURE_COORD_ARRAY_SIZE = GL_TEXTURE_COORD_ARRAY_SIZE,
// returns one value, the byte offset between consecutive elements in the texture coordinate array. The initial value is 0. See glTexCoordPointer.
TEXTURE_COORD_ARRAY_STRIDE = GL_TEXTURE_COORD_ARRAY_STRIDE,
// returns one value, the data type of the coordinates in the texture coordinate array. The initial value is GL_FLOAT. See glTexCoordPointer.
TEXTURE_COORD_ARRAY_TYPE = GL_TEXTURE_COORD_ARRAY_TYPE,
// returns a single boolean value indicating whether cube-mapped texture mapping is enabled. The initial value is GL_FALSE. See glTexImage2D.
TEXTURE_CUBE_MAP = GL_TEXTURE_CUBE_MAP,
// returns a single boolean value indicating whether automatic generation of the q texture coordinate is enabled. The initial value is GL_FALSE. See glTexGen.
TEXTURE_GEN_Q = GL_TEXTURE_GEN_Q,
// returns a single boolean value indicating whether automatic generation of the r texture coordinate is enabled. The initial value is GL_FALSE. See glTexGen.
TEXTURE_GEN_R = GL_TEXTURE_GEN_R,
// returns a single boolean value indicating whether automatic generation of the S texture coordinate is enabled. The initial value is GL_FALSE. See glTexGen.
TEXTURE_GEN_S = GL_TEXTURE_GEN_S,
// returns a single boolean value indicating whether automatic generation of the T texture coordinate is enabled. The initial value is GL_FALSE. See glTexGen.
TEXTURE_GEN_T = GL_TEXTURE_GEN_T,
// returns sixteen values: the texture matrix on the top of the texture matrix stack. Initially this matrix is the identity matrix. See glPushMatrix.
TEXTURE_MATRIX = GL_TEXTURE_MATRIX,
// returns one value, the number of matrices on the texture matrix stack. The initial value is 1. See glPushMatrix.
TEXTURE_STACK_DEPTH = GL_TEXTURE_STACK_DEPTH,
// returns 16 values, the elements of the color matrix in row-major order. See glLoadTransposeMatrix.
TRANSPOSE_COLOR_MATRIX = GL_TRANSPOSE_COLOR_MATRIX,
// returns 16 values, the elements of the modelview matrix in row-major order. See glLoadTransposeMatrix.
TRANSPOSE_MODELVIEW_MATRIX = GL_TRANSPOSE_MODELVIEW_MATRIX,
// returns 16 values, the elements of the projection matrix in row-major order. See glLoadTransposeMatrix.
TRANSPOSE_PROJECTION_MATRIX = GL_TRANSPOSE_PROJECTION_MATRIX,
// returns 16 values, the elements of the texture matrix in row-major order. See glLoadTransposeMatrix.
TRANSPOSE_TEXTURE_MATRIX = GL_TRANSPOSE_TEXTURE_MATRIX,
// returns one value, the byte alignment used for reading pixel data from memory. The initial value is 4. See glPixelStore.
UNPACK_ALIGNMENT = GL_UNPACK_ALIGNMENT,
// returns one value, the image height used for reading pixel data from memory. The initial is 0. See glPixelStore.
UNPACK_IMAGE_HEIGHT = GL_UNPACK_IMAGE_HEIGHT,
// returns a single boolean value indicating whether single-bit pixels being read from memory are read first from the least significant bit of each unsigned byte. The initial value is GL_FALSE. See glPixelStore.
UNPACK_LSB_FIRST = GL_UNPACK_LSB_FIRST,
// returns one value, the row length used for reading pixel data from memory. The initial value is 0. See glPixelStore.
UNPACK_ROW_LENGTH = GL_UNPACK_ROW_LENGTH,
// returns one value, the number of pixel images skipped before the first pixel is read from memory. The initial value is 0. See glPixelStore.
UNPACK_SKIP_IMAGES = GL_UNPACK_SKIP_IMAGES,
// returns one value, the number of pixel locations skipped before the first pixel is read from memory. The initial value is 0. See glPixelStore.
UNPACK_SKIP_PIXELS = GL_UNPACK_SKIP_PIXELS,
// returns one value, the number of rows of pixel locations skipped before the first pixel is read from memory. The initial value is 0. See glPixelStore.
UNPACK_SKIP_ROWS = GL_UNPACK_SKIP_ROWS,
// returns a single boolean value indicating whether the bytes of two-byte and four-byte pixel indices and components are swapped after being read from memory. The initial value is GL_FALSE. See glPixelStore.
UNPACK_SWAP_BYTES = GL_UNPACK_SWAP_BYTES,
// returns a single boolean value indicating whether the vertex array is enabled. The initial value is GL_FALSE. See glVertexPointer.
VERTEX_ARRAY = GL_VERTEX_ARRAY,
// returns a single value, the name of the buffer object associated with the vertex array. This buffer object would have been bound to the target GL_ARRAY_BUFFER at the time of the most recent call to glVertexPointer. If no buffer object was bound to this target, 0 is returned. The initial value is 0. See glBindBuffer.
VERTEX_ARRAY_BUFFER_BINDING = GL_VERTEX_ARRAY_BUFFER_BINDING,
// returns one value, the number of coordinates per vertex in the vertex array. The initial value is 4. See glVertexPointer.
VERTEX_ARRAY_SIZE = GL_VERTEX_ARRAY_SIZE,
// returns one value, the byte offset between consecutive vertices in the vertex array. The initial value is 0. See glVertexPointer.
VERTEX_ARRAY_STRIDE = GL_VERTEX_ARRAY_STRIDE,
// returns one value, the data type of each coordinate in the vertex array. The initial value is GL_FLOAT. See glVertexPointer.
VERTEX_ARRAY_TYPE = GL_VERTEX_ARRAY_TYPE,
// returns a single boolean value indicating whether vertex program point size mode is enabled. If enabled, and a vertex shader is active, then the point size is taken from the shader built-in gl_PointSize. If disabled, and a vertex shader is active, then the point size is taken from the point state as specified by glPointSize. The initial value is GL_FALSE.
VERTEX_PROGRAM_POINT_SIZE = GL_VERTEX_PROGRAM_POINT_SIZE,
// returns a single boolean value indicating whether vertex program two-sided color mode is enabled. If enabled, and a vertex shader is active, then the GL chooses the back color output for back-facing polygons, and the front color output for non-polygons and front-facing polygons. If disabled, and a vertex shader is active, then the front color output is always selected. The initial value is GL_FALSE.
VERTEX_PROGRAM_TWO_SIDE = GL_VERTEX_PROGRAM_TWO_SIDE,
// returns four values: the x and y window coordinates of the viewport, followed by its width and height. Initially the x and y window coordinates are both set to 0, and the width and height are set to the width and height of the window into which the GL will do its rendering. See glViewport.
VIEWPORT = GL_VIEWPORT,
// returns one value, the x pixel zoom factor. The initial value is 1. See glPixelZoom.
ZOOM_X = GL_ZOOM_X,
