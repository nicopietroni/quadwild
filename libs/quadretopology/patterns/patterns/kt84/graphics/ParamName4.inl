// data returns a single value indicating the active multitexture unit. The initial value is GL_TEXTURE0. See glActiveTexture.
ACTIVE_TEXTURE = GL_ACTIVE_TEXTURE,
// data returns a pair of values indicating the range of widths supported for aliased lines. See glLineWidth.
ALIASED_LINE_WIDTH_RANGE = GL_ALIASED_LINE_WIDTH_RANGE,
// returns a single value, the name of the buffer object currently bound to the target GL_ARRAY_BUFFER. If no buffer object is bound to this target, 0 is returned. The initial value is 0. See glBindBuffer.
ARRAY_BUFFER_BINDING = GL_ARRAY_BUFFER_BINDING,
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
// returns four values: the red, green, blue, and alpha values used to clear the color buffers. Integer values, if requested, are linearly mapped from the internal floating-point representation such that 1.0 returns the most positive representable integer value, and -1.0 returns the most negative representable integer value. The initial value is (0, 0, 0, 0). See glClearColor.
COLOR_CLEAR_VALUE = GL_COLOR_CLEAR_VALUE,
// returns a single boolean value indicating whether a fragment's RGBA color values are merged into the framebuffer using a logical operation. The initial value is GL_FALSE. See glLogicOp.
COLOR_LOGIC_OP = GL_COLOR_LOGIC_OP,
// returns four boolean values: the red, green, blue, and alpha write enables for the color buffers. The initial value is (GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE). See glColorMask.
COLOR_WRITEMASK = GL_COLOR_WRITEMASK,
// returns a list of symbolic constants of length GL_NUM_COMPRESSED_TEXTURE_FORMATS indicating which compressed texture formats are available. See glCompressedTexImage2D.
COMPRESSED_TEXTURE_FORMATS = GL_COMPRESSED_TEXTURE_FORMATS,
// returns one value, the maximum number of active shader storage blocks that may be accessed by a compute shader.
MAX_COMPUTE_SHADER_STORAGE_BLOCKS = GL_MAX_COMPUTE_SHADER_STORAGE_BLOCKS,
// returns one value, the maximum total number of active shader storage blocks that may be accessed by all active shaders.
MAX_COMBINED_SHADER_STORAGE_BLOCKS = GL_MAX_COMBINED_SHADER_STORAGE_BLOCKS,
// returns one value, the maximum number of uniform blocks per compute shader. The value must be at least 14. See glUniformBlockBinding.
MAX_COMPUTE_UNIFORM_BLOCKS = GL_MAX_COMPUTE_UNIFORM_BLOCKS,
// returns one value, the maximum supported texture image units that can be used to access texture maps from the compute shader. The value may be at least 16. See glActiveTexture.
MAX_COMPUTE_TEXTURE_IMAGE_UNITS = GL_MAX_COMPUTE_TEXTURE_IMAGE_UNITS,
// returns one value, the maximum number of individual floating-point, integer, or boolean values that can be held in uniform variable storage for a compute shader. The value must be at least 1024. See glUniform.
MAX_COMPUTE_UNIFORM_COMPONENTS = GL_MAX_COMPUTE_UNIFORM_COMPONENTS,
// returns a single value, the maximum number of atomic counters available to compute shaders.
MAX_COMPUTE_ATOMIC_COUNTERS = GL_MAX_COMPUTE_ATOMIC_COUNTERS,
// returns a single value, the maximum number of atomic counter buffers that may be accessed by a compute shader.
MAX_COMPUTE_ATOMIC_COUNTER_BUFFERS = GL_MAX_COMPUTE_ATOMIC_COUNTER_BUFFERS,
// returns one value, the number of words for compute shader uniform variables in all uniform blocks (including default). The value must be at least 1. See glUniform.
MAX_COMBINED_COMPUTE_UNIFORM_COMPONENTS = GL_MAX_COMBINED_COMPUTE_UNIFORM_COMPONENTS,
// returns one value, the number of invocations in a single local work group (i.e., the product of the three dimensions) that may be dispatched to a compute shader.
MAX_COMPUTE_WORK_GROUP_INVOCATIONS = GL_MAX_COMPUTE_WORK_GROUP_INVOCATIONS,
// Accepted by the indexed versions of glGet. the maximum number of work groups that may be dispatched to a compute shader. Indices 0, 1, and 2 correspond to the X, Y and Z dimensions, respectively.
MAX_COMPUTE_WORK_GROUP_COUNT = GL_MAX_COMPUTE_WORK_GROUP_COUNT,
// Accepted by the indexed versions of glGet. the maximum size of a work groups that may be used during compilation of a compute shader. Indices 0, 1, and 2 correspond to the X, Y and Z dimensions, respectively.
MAX_COMPUTE_WORK_GROUP_SIZE = GL_MAX_COMPUTE_WORK_GROUP_SIZE,
// returns a single value, the name of the buffer object currently bound to the target GL_DISPATCH_INDIRECT_BUFFER. If no buffer object is bound to this target, 0 is returned. The initial value is 0. See glBindBuffer.
DISPATCH_INDIRECT_BUFFER_BINDING = GL_DISPATCH_INDIRECT_BUFFER_BINDING,
// returns a single value, the maximum depth of the debug message group stack.
MAX_DEBUG_GROUP_STACK_DEPTH = GL_MAX_DEBUG_GROUP_STACK_DEPTH,
// returns a single value, the current depth of the debug message group stack.
DEBUG_GROUP_STACK_DEPTH = GL_DEBUG_GROUP_STACK_DEPTH,
// returns one value, the flags with which the context was created (such as debugging functionality).
CONTEXT_FLAGS = GL_CONTEXT_FLAGS,
// returns a single boolean value indicating whether polygon culling is enabled. The initial value is GL_FALSE. See glCullFace.
CULL_FACE = GL_CULL_FACE,
// returns one value, the name of the program object that is currently active, or 0 if no program object is active. See glUseProgram.
CURRENT_PROGRAM = GL_CURRENT_PROGRAM,
// returns one value, the value that is used to clear the depth buffer. Integer values, if requested, are linearly mapped from the internal floating-point representation such that 1.0 returns the most positive representable integer value, and -1.0 returns the most negative representable integer value. The initial value is 1. See glClearDepth.
DEPTH_CLEAR_VALUE = GL_DEPTH_CLEAR_VALUE,
// returns one value, the symbolic constant that indicates the depth comparison function. The initial value is GL_LESS. See glDepthFunc.
DEPTH_FUNC = GL_DEPTH_FUNC,
// returns two values: the near and far mapping limits for the depth buffer. Integer values, if requested, are linearly mapped from the internal floating-point representation such that 1.0 returns the most positive representable integer value, and -1.0 returns the most negative representable integer value. The initial value is (0, 1). See glDepthRange.
DEPTH_RANGE = GL_DEPTH_RANGE,
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
// returns one value, the name of the framebuffer object currently bound to the GL_DRAW_FRAMEBUFFER target. If the default framebuffer is bound, this value will be zero. The initial value is zero. See glBindFramebuffer.
DRAW_FRAMEBUFFER_BINDING = GL_DRAW_FRAMEBUFFER_BINDING,
// returns one value, the name of the framebuffer object currently bound to the GL_READ_FRAMEBUFFER target. If the default framebuffer is bound, this value will be zero. The initial value is zero. See glBindFramebuffer.
READ_FRAMEBUFFER_BINDING = GL_READ_FRAMEBUFFER_BINDING,
// returns a single value, the name of the buffer object currently bound to the target GL_ELEMENT_ARRAY_BUFFER. If no buffer object is bound to this target, 0 is returned. The initial value is 0. See glBindBuffer.
ELEMENT_ARRAY_BUFFER_BINDING = GL_ELEMENT_ARRAY_BUFFER_BINDING,
// returns one value, a symbolic constant indicating the mode of the derivative accuracy hint for fragment shaders. The initial value is GL_DONT_CARE. See glHint.
FRAGMENT_SHADER_DERIVATIVE_HINT = GL_FRAGMENT_SHADER_DERIVATIVE_HINT,
// returns a single GLenum value indicating the implementation's preferred pixel format. See glReadPixels.
IMPLEMENTATION_COLOR_READ_FORMAT = GL_IMPLEMENTATION_COLOR_READ_FORMAT,
// returns a single GLenum value indicating the implementation's preferred pixel type. See glReadPixels.
IMPLEMENTATION_COLOR_READ_TYPE = GL_IMPLEMENTATION_COLOR_READ_TYPE,
// returns a single boolean value indicating whether antialiasing of lines is enabled. The initial value is GL_FALSE. See glLineWidth.
LINE_SMOOTH = GL_LINE_SMOOTH,
// returns one value, a symbolic constant indicating the mode of the line antialiasing hint. The initial value is GL_DONT_CARE. See glHint.
LINE_SMOOTH_HINT = GL_LINE_SMOOTH_HINT,
// returns one value, the line width as specified with glLineWidth. The initial value is 1.
LINE_WIDTH = GL_LINE_WIDTH,
// returns one value, the implementation dependent specifc vertex of a primitive that is used to select the rendering layer. If the value returned is equivalent to GL_PROVOKING_VERTEX, then the vertex selection follows the convention specified by glProvokingVertex. If the value returned is equivalent to GL_FIRST_VERTEX_CONVENTION, then the selection is always taken from the first vertex in the primitive. If the value returned is equivalent to GL_LAST_VERTEX_CONVENTION, then the selection is always taken from the last vertex in the primitive. If the value returned is equivalent to GL_UNDEFINED_VERTEX, then the selection is not guaranteed to be taken from any specific vertex in the primitive.
LAYER_PROVOKING_VERTEX = GL_LAYER_PROVOKING_VERTEX,
// returns one value, a symbolic constant indicating the selected logic operation mode. The initial value is GL_COPY. See glLogicOp.
LOGIC_OP_MODE = GL_LOGIC_OP_MODE,
// returns one value, the major version number of the OpenGL API supported by the current context.
MAJOR_VERSION = GL_MAJOR_VERSION,
// returns one value, a rough estimate of the largest 3D texture that the GL can handle. The value must be at least 64. Use GL_PROXY_TEXTURE_3D to determine if a texture is too large. See glTexImage3D.
MAX_3D_TEXTURE_SIZE = GL_MAX_3D_TEXTURE_SIZE,
// returns one value. The value indicates the maximum number of layers allowed in an array texture, and must be at least 256. See glTexImage2D.
MAX_ARRAY_TEXTURE_LAYERS = GL_MAX_ARRAY_TEXTURE_LAYERS,
// returns one value, the maximum number of application-defined clipping distances. The value must be at least 8.
MAX_CLIP_DISTANCES = GL_MAX_CLIP_DISTANCES,
// returns one value, the maximum number of samples in a color multisample texture.
MAX_COLOR_TEXTURE_SAMPLES = GL_MAX_COLOR_TEXTURE_SAMPLES,
// returns a single value, the maximum number of atomic counters available to all active shaders.
MAX_COMBINED_ATOMIC_COUNTERS = GL_MAX_COMBINED_ATOMIC_COUNTERS,
// returns one value, the number of words for fragment shader uniform variables in all uniform blocks (including default). The value must be at least 1. See glUniform.
MAX_COMBINED_FRAGMENT_UNIFORM_COMPONENTS = GL_MAX_COMBINED_FRAGMENT_UNIFORM_COMPONENTS,
// returns one value, the number of words for geometry shader uniform variables in all uniform blocks (including default). The value must be at least 1. See glUniform.
MAX_COMBINED_GEOMETRY_UNIFORM_COMPONENTS = GL_MAX_COMBINED_GEOMETRY_UNIFORM_COMPONENTS,
// returns one value, the maximum supported texture image units that can be used to access texture maps from the vertex shader and the fragment processor combined. If both the vertex shader and the fragment processing stage access the same texture image unit, then that counts as using two texture image units against this limit. The value must be at least 48. See glActiveTexture.
MAX_COMBINED_TEXTURE_IMAGE_UNITS = GL_MAX_COMBINED_TEXTURE_IMAGE_UNITS,
// returns one value, the maximum number of uniform blocks per program. The value must be at least 70. See glUniformBlockBinding.
MAX_COMBINED_UNIFORM_BLOCKS = GL_MAX_COMBINED_UNIFORM_BLOCKS,
// returns one value, the number of words for vertex shader uniform variables in all uniform blocks (including default). The value must be at least 1. See glUniform.
MAX_COMBINED_VERTEX_UNIFORM_COMPONENTS = GL_MAX_COMBINED_VERTEX_UNIFORM_COMPONENTS,
// returns one value. The value gives a rough estimate of the largest cube-map texture that the GL can handle. The value must be at least 1024. Use GL_PROXY_TEXTURE_CUBE_MAP to determine if a texture is too large. See glTexImage2D.
MAX_CUBE_MAP_TEXTURE_SIZE = GL_MAX_CUBE_MAP_TEXTURE_SIZE,
// returns one value, the maximum number of samples in a multisample depth or depth-stencil texture.
MAX_DEPTH_TEXTURE_SAMPLES = GL_MAX_DEPTH_TEXTURE_SAMPLES,
// returns one value, the maximum number of simultaneous outputs that may be written in a fragment shader. The value must be at least 8. See glDrawBuffers.
MAX_DRAW_BUFFERS = GL_MAX_DRAW_BUFFERS,
// returns one value, the maximum number of active draw buffers when using dual-source blending. The value must be at least 1. See glBlendFunc and glBlendFuncSeparate.
MAX_DUAL_SOURCE_DRAW_BUFFERS = GL_MAX_DUAL_SOURCE_DRAW_BUFFERS,
// returns one value, the recommended maximum number of vertex array indices. See glDrawRangeElements.
MAX_ELEMENTS_INDICES = GL_MAX_ELEMENTS_INDICES,
// returns one value, the recommended maximum number of vertex array vertices. See glDrawRangeElements.
MAX_ELEMENTS_VERTICES = GL_MAX_ELEMENTS_VERTICES,
// returns a single value, the maximum number of atomic counters available to fragment shaders.
MAX_FRAGMENT_ATOMIC_COUNTERS = GL_MAX_FRAGMENT_ATOMIC_COUNTERS,
// returns one value, the maximum number of active shader storage blocks that may be accessed by a fragment shader.
MAX_FRAGMENT_SHADER_STORAGE_BLOCKS = GL_MAX_FRAGMENT_SHADER_STORAGE_BLOCKS,
// returns one value, the maximum number of components of the inputs read by the fragment shader, which must be at least 128.
MAX_FRAGMENT_INPUT_COMPONENTS = GL_MAX_FRAGMENT_INPUT_COMPONENTS,
// returns one value, the maximum number of individual floating-point, integer, or boolean values that can be held in uniform variable storage for a fragment shader. The value must be at least 1024. See glUniform.
MAX_FRAGMENT_UNIFORM_COMPONENTS = GL_MAX_FRAGMENT_UNIFORM_COMPONENTS,
// returns one value, the maximum number of individual 4-vectors of floating-point, integer, or boolean values that can be held in uniform variable storage for a fragment shader. The value is equal to the value of GL_MAX_FRAGMENT_UNIFORM_COMPONENTS divided by 4 and must be at least 256. See glUniform.
MAX_FRAGMENT_UNIFORM_VECTORS = GL_MAX_FRAGMENT_UNIFORM_VECTORS,
// returns one value, the maximum number of uniform blocks per fragment shader. The value must be at least 12. See glUniformBlockBinding.
MAX_FRAGMENT_UNIFORM_BLOCKS = GL_MAX_FRAGMENT_UNIFORM_BLOCKS,
// returns one value, the maximum width for a framebuffer that has no attachments, which must be at least 16384. See glFramebufferParameter.
MAX_FRAMEBUFFER_WIDTH = GL_MAX_FRAMEBUFFER_WIDTH,
// returns one value, the maximum height for a framebuffer that has no attachments, which must be at least 16384. See glFramebufferParameter.
MAX_FRAMEBUFFER_HEIGHT = GL_MAX_FRAMEBUFFER_HEIGHT,
// returns one value, the maximum number of layers for a framebuffer that has no attachments, which must be at least 2048. See glFramebufferParameter.
MAX_FRAMEBUFFER_LAYERS = GL_MAX_FRAMEBUFFER_LAYERS,
// returns one value, the maximum samples in a framebuffer that has no attachments, which must be at least 4. See glFramebufferParameter.
MAX_FRAMEBUFFER_SAMPLES = GL_MAX_FRAMEBUFFER_SAMPLES,
// returns a single value, the maximum number of atomic counters available to geometry shaders.
MAX_GEOMETRY_ATOMIC_COUNTERS = GL_MAX_GEOMETRY_ATOMIC_COUNTERS,
// returns one value, the maximum number of active shader storage blocks that may be accessed by a geometry shader.
MAX_GEOMETRY_SHADER_STORAGE_BLOCKS = GL_MAX_GEOMETRY_SHADER_STORAGE_BLOCKS,
// returns one value, the maximum number of components of inputs read by a geometry shader, which must be at least 64.
MAX_GEOMETRY_INPUT_COMPONENTS = GL_MAX_GEOMETRY_INPUT_COMPONENTS,
// returns one value, the maximum number of components of outputs written by a geometry shader, which must be at least 128.
MAX_GEOMETRY_OUTPUT_COMPONENTS = GL_MAX_GEOMETRY_OUTPUT_COMPONENTS,
// returns one value, the maximum supported texture image units that can be used to access texture maps from the geometry shader. The value must be at least 16. See glActiveTexture.
MAX_GEOMETRY_TEXTURE_IMAGE_UNITS = GL_MAX_GEOMETRY_TEXTURE_IMAGE_UNITS,
// returns one value, the maximum number of uniform blocks per geometry shader. The value must be at least 12. See glUniformBlockBinding.
MAX_GEOMETRY_UNIFORM_BLOCKS = GL_MAX_GEOMETRY_UNIFORM_BLOCKS,
// returns one value, the maximum number of individual floating-point, integer, or boolean values that can be held in uniform variable storage for a geometry shader. The value must be at least 1024. See glUniform.
MAX_GEOMETRY_UNIFORM_COMPONENTS = GL_MAX_GEOMETRY_UNIFORM_COMPONENTS,
// returns one value, the maximum number of samples supported in integer format multisample buffers.
MAX_INTEGER_SAMPLES = GL_MAX_INTEGER_SAMPLES,
// returns one value, the minimum alignment in basic machine units of pointers returned fromglMapBuffer and glMapBufferRange. This value must be a power of two and must be at least 64.
MIN_MAP_BUFFER_ALIGNMENT = GL_MIN_MAP_BUFFER_ALIGNMENT,
// returns one value, the maximum length of a label that may be assigned to an object. See glObjectLabel and glObjectPtrLabel.
MAX_LABEL_LENGTH = GL_MAX_LABEL_LENGTH,
// returns one value, the maximum texel offset allowed in a texture lookup, which must be at least 7.
MAX_PROGRAM_TEXEL_OFFSET = GL_MAX_PROGRAM_TEXEL_OFFSET,
// returns one value, the minimum texel offset allowed in a texture lookup, which must be at most -8.
MIN_PROGRAM_TEXEL_OFFSET = GL_MIN_PROGRAM_TEXEL_OFFSET,
// returns one value. The value gives a rough estimate of the largest rectangular texture that the GL can handle. The value must be at least 1024. Use GL_PROXY_TEXTURE_RECTANGLE to determine if a texture is too large. See glTexImage2D.
MAX_RECTANGLE_TEXTURE_SIZE = GL_MAX_RECTANGLE_TEXTURE_SIZE,
// returns one value. The value indicates the maximum supported size for renderbuffers. See glFramebufferRenderbuffer.
MAX_RENDERBUFFER_SIZE = GL_MAX_RENDERBUFFER_SIZE,
// returns one value, the maximum number of sample mask words.
MAX_SAMPLE_MASK_WORDS = GL_MAX_SAMPLE_MASK_WORDS,
// returns one value, the maximum glWaitSync timeout interval.
MAX_SERVER_WAIT_TIMEOUT = GL_MAX_SERVER_WAIT_TIMEOUT,
// returns one value, the maximum number of shader storage buffer binding points on the context, which must be at least 8.
MAX_SHADER_STORAGE_BUFFER_BINDINGS = GL_MAX_SHADER_STORAGE_BUFFER_BINDINGS,
// returns a single value, the maximum number of atomic counters available to tessellation control shaders.
MAX_TESS_CONTROL_ATOMIC_COUNTERS = GL_MAX_TESS_CONTROL_ATOMIC_COUNTERS,
// returns a single value, the maximum number of atomic counters available to tessellation evaluation shaders.
MAX_TESS_EVALUATION_ATOMIC_COUNTERS = GL_MAX_TESS_EVALUATION_ATOMIC_COUNTERS,
// returns one value, the maximum number of active shader storage blocks that may be accessed by a tessellation control shader.
MAX_TESS_CONTROL_SHADER_STORAGE_BLOCKS = GL_MAX_TESS_CONTROL_SHADER_STORAGE_BLOCKS,
// returns one value, the maximum number of active shader storage blocks that may be accessed by a tessellation evaluation shader.
MAX_TESS_EVALUATION_SHADER_STORAGE_BLOCKS = GL_MAX_TESS_EVALUATION_SHADER_STORAGE_BLOCKS,
// returns one value. The value gives the maximum number of texels allowed in the texel array of a texture buffer object. Value must be at least 65536.
MAX_TEXTURE_BUFFER_SIZE = GL_MAX_TEXTURE_BUFFER_SIZE,
// returns one value, the maximum supported texture image units that can be used to access texture maps from the fragment shader. The value must be at least 16. See glActiveTexture.
MAX_TEXTURE_IMAGE_UNITS = GL_MAX_TEXTURE_IMAGE_UNITS,
// returns one value, the maximum, absolute value of the texture level-of-detail bias. The value must be at least 2.0.
MAX_TEXTURE_LOD_BIAS = GL_MAX_TEXTURE_LOD_BIAS,
// returns one value. The value gives a rough estimate of the largest texture that the GL can handle. The value must be at least 1024. Use a proxy texture target such as GL_PROXY_TEXTURE_1D or GL_PROXY_TEXTURE_2D to determine if a texture is too large. See glTexImage1D and glTexImage2D.
MAX_TEXTURE_SIZE = GL_MAX_TEXTURE_SIZE,
// returns one value, the maximum number of uniform buffer binding points on the context, which must be at least 36.
MAX_UNIFORM_BUFFER_BINDINGS = GL_MAX_UNIFORM_BUFFER_BINDINGS,
// returns one value, the maximum size in basic machine units of a uniform block, which must be at least 16384.
MAX_UNIFORM_BLOCK_SIZE = GL_MAX_UNIFORM_BLOCK_SIZE,
// returns one value, the maximum number of explicitly assignable uniform locations, which must be at least 1024.
MAX_UNIFORM_LOCATIONS = GL_MAX_UNIFORM_LOCATIONS,
// returns one value, the number components for varying variables, which must be at least 60.
MAX_VARYING_COMPONENTS = GL_MAX_VARYING_COMPONENTS,
// returns one value, the number 4-vectors for varying variables, which is equal to the value of GL_MAX_VARYING_COMPONENTS and must be at least 15.
MAX_VARYING_VECTORS = GL_MAX_VARYING_VECTORS,
// returns one value, the maximum number of interpolators available for processing varying variables used by vertex and fragment shaders. This value represents the number of individual floating-point values that can be interpolated; varying variables declared as vectors, matrices, and arrays will all consume multiple interpolators. The value must be at least 32.
MAX_VARYING_FLOATS = GL_MAX_VARYING_FLOATS,
// returns a single value, the maximum number of atomic counters available to vertex shaders.
MAX_VERTEX_ATOMIC_COUNTERS = GL_MAX_VERTEX_ATOMIC_COUNTERS,
// returns one value, the maximum number of 4-component generic vertex attributes accessible to a vertex shader. The value must be at least 16. See glVertexAttrib.
MAX_VERTEX_ATTRIBS = GL_MAX_VERTEX_ATTRIBS,
// returns one value, the maximum number of active shader storage blocks that may be accessed by a vertex shader.
MAX_VERTEX_SHADER_STORAGE_BLOCKS = GL_MAX_VERTEX_SHADER_STORAGE_BLOCKS,
// returns one value, the maximum supported texture image units that can be used to access texture maps from the vertex shader. The value may be at least 16. See glActiveTexture.
MAX_VERTEX_TEXTURE_IMAGE_UNITS = GL_MAX_VERTEX_TEXTURE_IMAGE_UNITS,
// returns one value, the maximum number of individual floating-point, integer, or boolean values that can be held in uniform variable storage for a vertex shader. The value must be at least 1024. See glUniform.
MAX_VERTEX_UNIFORM_COMPONENTS = GL_MAX_VERTEX_UNIFORM_COMPONENTS,
// returns one value, the maximum number of 4-vectors that may be held in uniform variable storage for the vertex shader. The value of GL_MAX_VERTEX_UNIFORM_VECTORS is equal to the value of GL_MAX_VERTEX_UNIFORM_COMPONENTS and must be at least 256.
MAX_VERTEX_UNIFORM_VECTORS = GL_MAX_VERTEX_UNIFORM_VECTORS,
// returns one value, the maximum number of components of output written by a vertex shader, which must be at least 64.
MAX_VERTEX_OUTPUT_COMPONENTS = GL_MAX_VERTEX_OUTPUT_COMPONENTS,
// returns one value, the maximum number of uniform blocks per vertex shader. The value must be at least 12. See glUniformBlockBinding.
MAX_VERTEX_UNIFORM_BLOCKS = GL_MAX_VERTEX_UNIFORM_BLOCKS,
// returns two values: the maximum supported width and height of the viewport. These must be at least as large as the visible dimensions of the display being rendered to. See glViewport.
MAX_VIEWPORT_DIMS = GL_MAX_VIEWPORT_DIMS,
// returns one value, the maximum number of simultaneous viewports that are supported. The value must be at least 16. See glViewportIndexed.
MAX_VIEWPORTS = GL_MAX_VIEWPORTS,
// returns one value, the minor version number of the OpenGL API supported by the current context.
MINOR_VERSION = GL_MINOR_VERSION,
// returns a single integer value indicating the number of available compressed texture formats. The minimum value is 4. See glCompressedTexImage2D.
NUM_COMPRESSED_TEXTURE_FORMATS = GL_NUM_COMPRESSED_TEXTURE_FORMATS,
// returns one value, the number of extensions supported by the GL implementation for the current context. See glGetString.
NUM_EXTENSIONS = GL_NUM_EXTENSIONS,
// returns one value, the number of program binary formats supported by the implementation.
NUM_PROGRAM_BINARY_FORMATS = GL_NUM_PROGRAM_BINARY_FORMATS,
// returns one value, the number of binary shader formats supported by the implementation. If this value is greater than zero, then the implementation supports loading binary shaders. If it is zero, then the loading of binary shaders by the implementation is not supported.
NUM_SHADER_BINARY_FORMATS = GL_NUM_SHADER_BINARY_FORMATS,
// returns one value, the byte alignment used for writing pixel to memory. The initial value is 4. See glPixelStore.
PACK_ALIGNMENT = GL_PACK_ALIGNMENT,
// returns one value, the image height used for writing pixel to memory. The initial value is 0. See glPixelStore.
PACK_IMAGE_HEIGHT = GL_PACK_IMAGE_HEIGHT,
// returns a single boolean value indicating whether single-bit pixels being written to memory are written first to the least significant bit of each unsigned byte. The initial value is GL_FALSE. See glPixelStore.
PACK_LSB_FIRST = GL_PACK_LSB_FIRST,
// returns one value, the row length used for writing pixel to memory. The initial value is 0. See glPixelStore.
PACK_ROW_LENGTH = GL_PACK_ROW_LENGTH,
// returns one value, the number of pixel images skipped before the first pixel is written into memory. The initial value is 0. See glPixelStore.
PACK_SKIP_IMAGES = GL_PACK_SKIP_IMAGES,
// returns one value, the number of pixel locations skipped before the first pixel is written into memory. The initial value is 0. See glPixelStore.
PACK_SKIP_PIXELS = GL_PACK_SKIP_PIXELS,
// returns one value, the number of rows of pixel locations skipped before the first pixel is written into memory. The initial value is 0. See glPixelStore.
PACK_SKIP_ROWS = GL_PACK_SKIP_ROWS,
// returns a single boolean value indicating whether the bytes of two-byte and four-byte pixel indices and components are swapped before being written to memory. The initial value is GL_FALSE. See glPixelStore.
PACK_SWAP_BYTES = GL_PACK_SWAP_BYTES,
// returns a single value, the name of the buffer object currently bound to the target GL_PIXEL_PACK_BUFFER. If no buffer object is bound to this target, 0 is returned. The initial value is 0. See glBindBuffer.
PIXEL_PACK_BUFFER_BINDING = GL_PIXEL_PACK_BUFFER_BINDING,
// returns a single value, the name of the buffer object currently bound to the target GL_PIXEL_UNPACK_BUFFER. If no buffer object is bound to this target, 0 is returned. The initial value is 0. See glBindBuffer.
PIXEL_UNPACK_BUFFER_BINDING = GL_PIXEL_UNPACK_BUFFER_BINDING,
// returns one value, the point size threshold for determining the point size. See glPointParameter.
POINT_FADE_THRESHOLD_SIZE = GL_POINT_FADE_THRESHOLD_SIZE,
// returns one value, the current primitive restart index. The initial value is 0. See glPrimitiveRestartIndex.
PRIMITIVE_RESTART_INDEX = GL_PRIMITIVE_RESTART_INDEX,
// an array of GL_NUM_PROGRAM_BINARY_FORMATS values, indicating the proram binary formats supported by the implementation.
PROGRAM_BINARY_FORMATS = GL_PROGRAM_BINARY_FORMATS,
// a single value, the name of the currently bound program pipeline object, or zero if no program pipeline object is bound. See glBindProgramPipeline.
PROGRAM_PIPELINE_BINDING = GL_PROGRAM_PIPELINE_BINDING,
// returns a single boolean value indicating whether vertex program point size mode is enabled. If enabled, then the point size is taken from the shader built-in gl_PointSize. If disabled, then the point size is taken from the point state as specified by glPointSize. The initial value is GL_FALSE.
PROGRAM_POINT_SIZE = GL_PROGRAM_POINT_SIZE,
// returns one value, the currently selected provoking vertex convention. The initial value is GL_LAST_VERTEX_CONVENTION. See glProvokingVertex.
PROVOKING_VERTEX = GL_PROVOKING_VERTEX,
// returns one value, the point size as specified by glPointSize. The initial value is 1.
POINT_SIZE = GL_POINT_SIZE,
// returns one value, the size difference between adjacent supported sizes for antialiased points. See glPointSize.
POINT_SIZE_GRANULARITY = GL_POINT_SIZE_GRANULARITY,
// returns two values: the smallest and largest supported sizes for antialiased points. The smallest size must be at most 1, and the largest size must be at least 1. See glPointSize.
POINT_SIZE_RANGE = GL_POINT_SIZE_RANGE,
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
// returns one value, a symbolic constant indicating which color buffer is selected for reading. The initial value is GL_BACK if there is a back buffer, otherwise it is GL_FRONT. See glReadPixels.
READ_BUFFER = GL_READ_BUFFER,
// returns a single value, the name of the renderbuffer object currently bound to the target GL_RENDERBUFFER. If no renderbuffer object is bound to this target, 0 is returned. The initial value is 0. See glBindRenderbuffer.
RENDERBUFFER_BINDING = GL_RENDERBUFFER_BINDING,
// returns a single integer value indicating the number of sample buffers associated with the framebuffer. See glSampleCoverage.
SAMPLE_BUFFERS = GL_SAMPLE_BUFFERS,
// returns a single positive floating-point value indicating the current sample coverage value. See glSampleCoverage.
SAMPLE_COVERAGE_VALUE = GL_SAMPLE_COVERAGE_VALUE,
// returns a single boolean value indicating if the temporary coverage value should be inverted. See glSampleCoverage.
SAMPLE_COVERAGE_INVERT = GL_SAMPLE_COVERAGE_INVERT,
// returns a single value, the name of the sampler object currently bound to the active texture unit. The initial value is 0. See glBindSampler.
SAMPLER_BINDING = GL_SAMPLER_BINDING,
// returns a single integer value indicating the coverage mask size. See glSampleCoverage.
SAMPLES = GL_SAMPLES,
// returns four values: the x and y window coordinates of the scissor box, followed by its width and height. Initially the x and y window coordinates are both 0 and the width and height are set to the size of the window. See glScissor.
SCISSOR_BOX = GL_SCISSOR_BOX,
// returns a single boolean value indicating whether scissoring is enabled. The initial value is GL_FALSE. See glScissor.
SCISSOR_TEST = GL_SCISSOR_TEST,
// returns a single boolean value indicating whether an online shader compiler is present in the implementation. All desktop OpenGL implementations must support online shader compilations, and therefore the value of GL_SHADER_COMPILER will always be GL_TRUE.
SHADER_COMPILER = GL_SHADER_COMPILER,
// When used with non-indexed variants of glGet (such as glGetIntegerv), returns a single value, the name of the buffer object currently bound to the target GL_SHADER_STORAGE_BUFFER. If no buffer object is bound to this target, 0 is returned. When used with indexed variants of glGet (such as glGetIntegeri_v), returns a single value, the name of the buffer object bound to the indexed shader storage buffer binding points. The initial value is 0 for all targets. See glBindBuffer, glBindBufferBase, and glBindBufferRange.
SHADER_STORAGE_BUFFER_BINDING = GL_SHADER_STORAGE_BUFFER_BINDING,
// returns a single value, the minimum required alignment for shader storage buffer sizes and offset. The initial value is 1. See glShaderStorateBlockBinding.
SHADER_STORAGE_BUFFER_OFFSET_ALIGNMENT = GL_SHADER_STORAGE_BUFFER_OFFSET_ALIGNMENT,
// When used with indexed variants of glGet (such as glGetInteger64i_v), returns a single value, the start offset of the binding range for each indexed shader storage buffer binding. The initial value is 0 for all bindings. See glBindBufferRange.
SHADER_STORAGE_BUFFER_START = GL_SHADER_STORAGE_BUFFER_START,
// When used with indexed variants of glGet (such as glGetInteger64i_v), returns a single value, the size of the binding range for each indexed shader storage buffer binding. The initial value is 0 for all bindings. See glBindBufferRange.
SHADER_STORAGE_BUFFER_SIZE = GL_SHADER_STORAGE_BUFFER_SIZE,
// returns a pair of values indicating the range of widths supported for smooth (antialiased) lines. See glLineWidth.
SMOOTH_LINE_WIDTH_RANGE = GL_SMOOTH_LINE_WIDTH_RANGE,
// returns a single value indicating the level of quantization applied to smooth line width parameters.
SMOOTH_LINE_WIDTH_GRANULARITY = GL_SMOOTH_LINE_WIDTH_GRANULARITY,
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
// returns one value, the index to which the stencil bitplanes are cleared. The initial value is 0. See glClearStencil.
STENCIL_CLEAR_VALUE = GL_STENCIL_CLEAR_VALUE,
// returns one value, a symbolic constant indicating what action is taken when the stencil test fails. The initial value is GL_KEEP. See glStencilOp. This stencil state only affects non-polygons and front-facing polygons. Back-facing polygons use separate stencil state. See glStencilOpSeparate.
STENCIL_FAIL = GL_STENCIL_FAIL,
// returns one value, a symbolic constant indicating what function is used to compare the stencil reference value with the stencil buffer value. The initial value is GL_ALWAYS. See glStencilFunc. This stencil state only affects non-polygons and front-facing polygons. Back-facing polygons use separate stencil state. See glStencilFuncSeparate.
STENCIL_FUNC = GL_STENCIL_FUNC,
// returns one value, a symbolic constant indicating what action is taken when the stencil test passes, but the depth test fails. The initial value is GL_KEEP. See glStencilOp. This stencil state only affects non-polygons and front-facing polygons. Back-facing polygons use separate stencil state. See glStencilOpSeparate.
STENCIL_PASS_DEPTH_FAIL = GL_STENCIL_PASS_DEPTH_FAIL,
// returns one value, a symbolic constant indicating what action is taken when the stencil test passes and the depth test passes. The initial value is GL_KEEP. See glStencilOp. This stencil state only affects non-polygons and front-facing polygons. Back-facing polygons use separate stencil state. See glStencilOpSeparate.
STENCIL_PASS_DEPTH_PASS = GL_STENCIL_PASS_DEPTH_PASS,
// returns one value, the reference value that is compared with the contents of the stencil buffer. The initial value is 0. See glStencilFunc. This stencil state only affects non-polygons and front-facing polygons. Back-facing polygons use separate stencil state. See glStencilFuncSeparate.
STENCIL_REF = GL_STENCIL_REF,
// returns a single boolean value indicating whether stencil testing of fragments is enabled. The initial value is GL_FALSE. See glStencilFunc and glStencilOp.
STENCIL_TEST = GL_STENCIL_TEST,
// returns one value, the mask that is used to mask both the stencil reference value and the stencil buffer value before they are compared. The initial value is all 1's. See glStencilFunc. This stencil state only affects non-polygons and front-facing polygons. Back-facing polygons use separate stencil state. See glStencilFuncSeparate.
STENCIL_VALUE_MASK = GL_STENCIL_VALUE_MASK,
// returns one value, the mask that controls writing of the stencil bitplanes. The initial value is all 1's. See glStencilMask. This stencil state only affects non-polygons and front-facing polygons. Back-facing polygons use separate stencil state. See glStencilMaskSeparate.
STENCIL_WRITEMASK = GL_STENCIL_WRITEMASK,
// returns a single boolean value indicating whether stereo buffers (left and right) are supported.
STEREO = GL_STEREO,
// returns one value, an estimate of the number of bits of subpixel resolution that are used to position rasterized geometry in window coordinates. The value must be at least 4.
SUBPIXEL_BITS = GL_SUBPIXEL_BITS,
// returns a single value, the name of the texture currently bound to the target GL_TEXTURE_1D. The initial value is 0. See glBindTexture.
TEXTURE_BINDING_1D = GL_TEXTURE_BINDING_1D,
// returns a single value, the name of the texture currently bound to the target GL_TEXTURE_1D_ARRAY. The initial value is 0. See glBindTexture.
TEXTURE_BINDING_1D_ARRAY = GL_TEXTURE_BINDING_1D_ARRAY,
// returns a single value, the name of the texture currently bound to the target GL_TEXTURE_2D. The initial value is 0. See glBindTexture.
TEXTURE_BINDING_2D = GL_TEXTURE_BINDING_2D,
// returns a single value, the name of the texture currently bound to the target GL_TEXTURE_2D_ARRAY. The initial value is 0. See glBindTexture.
TEXTURE_BINDING_2D_ARRAY = GL_TEXTURE_BINDING_2D_ARRAY,
// returns a single value, the name of the texture currently bound to the target GL_TEXTURE_2D_MULTISAMPLE. The initial value is 0. See glBindTexture.
TEXTURE_BINDING_2D_MULTISAMPLE = GL_TEXTURE_BINDING_2D_MULTISAMPLE,
// returns a single value, the name of the texture currently bound to the target GL_TEXTURE_2D_MULTISAMPLE_ARRAY. The initial value is 0. See glBindTexture.
TEXTURE_BINDING_2D_MULTISAMPLE_ARRAY = GL_TEXTURE_BINDING_2D_MULTISAMPLE_ARRAY,
// returns a single value, the name of the texture currently bound to the target GL_TEXTURE_3D. The initial value is 0. See glBindTexture.
TEXTURE_BINDING_3D = GL_TEXTURE_BINDING_3D,
// returns a single value, the name of the texture currently bound to the target GL_TEXTURE_BUFFER. The initial value is 0. See glBindTexture.
TEXTURE_BINDING_BUFFER = GL_TEXTURE_BINDING_BUFFER,
// returns a single value, the name of the texture currently bound to the target GL_TEXTURE_CUBE_MAP. The initial value is 0. See glBindTexture.
TEXTURE_BINDING_CUBE_MAP = GL_TEXTURE_BINDING_CUBE_MAP,
// returns a single value, the name of the texture currently bound to the target GL_TEXTURE_RECTANGLE. The initial value is 0. See glBindTexture.
TEXTURE_BINDING_RECTANGLE = GL_TEXTURE_BINDING_RECTANGLE,
// returns a single value indicating the mode of the texture compression hint. The initial value is GL_DONT_CARE.
TEXTURE_COMPRESSION_HINT = GL_TEXTURE_COMPRESSION_HINT,
// returns a single value, the name of the buffer object currently bound to the GL_TEXTURE_BUFFER buffer binding point. The initial value is 0. See glBindBuffer.
TEXTURE_BINDING_BUFFER = GL_TEXTURE_BINDING_BUFFER,
// returns a single value, the minimum required alignment for texture buffer sizes and offset. The initial value is 1. See glUniformBlockBinding.
TEXTURE_BUFFER_OFFSET_ALIGNMENT = GL_TEXTURE_BUFFER_OFFSET_ALIGNMENT,
// returns a single value, the 64-bit value of the current GL time. See glQueryCounter.
TIMESTAMP = GL_TIMESTAMP,
// When used with non-indexed variants of glGet (such as glGetIntegerv), returns a single value, the name of the buffer object currently bound to the target GL_TRANSFORM_FEEDBACK_BUFFER. If no buffer object is bound to this target, 0 is returned. When used with indexed variants of glGet (such as glGetIntegeri_v), returns a single value, the name of the buffer object bound to the indexed transform feedback attribute stream. The initial value is 0 for all targets. See glBindBuffer, glBindBufferBase, and glBindBufferRange.
TRANSFORM_FEEDBACK_BUFFER_BINDING = GL_TRANSFORM_FEEDBACK_BUFFER_BINDING,
// When used with indexed variants of glGet (such as glGetInteger64i_v), returns a single value, the start offset of the binding range for each transform feedback attribute stream. The initial value is 0 for all streams. See glBindBufferRange.
TRANSFORM_FEEDBACK_BUFFER_START = GL_TRANSFORM_FEEDBACK_BUFFER_START,
// When used with indexed variants of glGet (such as glGetInteger64i_v), returns a single value, the size of the binding range for each transform feedback attribute stream. The initial value is 0 for all streams. See glBindBufferRange.
TRANSFORM_FEEDBACK_BUFFER_SIZE = GL_TRANSFORM_FEEDBACK_BUFFER_SIZE,
// When used with non-indexed variants of glGet (such as glGetIntegerv), returns a single value, the name of the buffer object currently bound to the target GL_UNIFORM_BUFFER. If no buffer object is bound to this target, 0 is returned. When used with indexed variants of glGet (such as glGetIntegeri_v), returns a single value, the name of the buffer object bound to the indexed uniform buffer binding point. The initial value is 0 for all targets. See glBindBuffer, glBindBufferBase, and glBindBufferRange.
UNIFORM_BUFFER_BINDING = GL_UNIFORM_BUFFER_BINDING,
// returns a single value, the minimum required alignment for uniform buffer sizes and offset. The initial value is 1. See glUniformBlockBinding.
UNIFORM_BUFFER_OFFSET_ALIGNMENT = GL_UNIFORM_BUFFER_OFFSET_ALIGNMENT,
// When used with indexed variants of glGet (such as glGetInteger64i_v), returns a single value, the size of the binding range for each indexed uniform buffer binding. The initial value is 0 for all bindings. See glBindBufferRange.
UNIFORM_BUFFER_SIZE = GL_UNIFORM_BUFFER_SIZE,
// When used with indexed variants of glGet (such as glGetInteger64i_v), returns a single value, the start offset of the binding range for each indexed uniform buffer binding. The initial value is 0 for all bindings. See glBindBufferRange.
UNIFORM_BUFFER_START = GL_UNIFORM_BUFFER_START,
// returns one value, the byte alignment used for reading pixel from memory. The initial value is 4. See glPixelStore.
UNPACK_ALIGNMENT = GL_UNPACK_ALIGNMENT,
// returns one value, the image height used for reading pixel from memory. The initial is 0. See glPixelStore.
UNPACK_IMAGE_HEIGHT = GL_UNPACK_IMAGE_HEIGHT,
// returns a single boolean value indicating whether single-bit pixels being read from memory are read first from the least significant bit of each unsigned byte. The initial value is GL_FALSE. See glPixelStore.
UNPACK_LSB_FIRST = GL_UNPACK_LSB_FIRST,
// returns one value, the row length used for reading pixel from memory. The initial value is 0. See glPixelStore.
UNPACK_ROW_LENGTH = GL_UNPACK_ROW_LENGTH,
// returns one value, the number of pixel images skipped before the first pixel is read from memory. The initial value is 0. See glPixelStore.
UNPACK_SKIP_IMAGES = GL_UNPACK_SKIP_IMAGES,
// returns one value, the number of pixel locations skipped before the first pixel is read from memory. The initial value is 0. See glPixelStore.
UNPACK_SKIP_PIXELS = GL_UNPACK_SKIP_PIXELS,
// returns one value, the number of rows of pixel locations skipped before the first pixel is read from memory. The initial value is 0. See glPixelStore.
UNPACK_SKIP_ROWS = GL_UNPACK_SKIP_ROWS,
// returns a single boolean value indicating whether the bytes of two-byte and four-byte pixel indices and components are swapped after being read from memory. The initial value is GL_FALSE. See glPixelStore.
UNPACK_SWAP_BYTES = GL_UNPACK_SWAP_BYTES,
// returns a single value, the name of the vertex array object currently bound to the context. If no vertex array object is bound to the context, 0 is returned. The initial value is 0. See glBindVertexArray.
VERTEX_ARRAY_BINDING = GL_VERTEX_ARRAY_BINDING,
// Accepted by the indexed forms. returns a single integer value representing the instance step divisor of the first element in the bound buffer's store for vertex attribute bound to index.
VERTEX_BINDING_DIVISOR = GL_VERTEX_BINDING_DIVISOR,
// Accepted by the indexed forms. returns a single integer value representing the byte offset of the first element in the bound buffer's store for vertex attribute bound to index.
VERTEX_BINDING_OFFSET = GL_VERTEX_BINDING_OFFSET,
// Accepted by the indexed forms. returns a single integer value representing the byte offset between the start of each element in the bound buffer's store for vertex attribute bound to index.
VERTEX_BINDING_STRIDE = GL_VERTEX_BINDING_STRIDE,
// returns a single integer value containing the maximum offset that may be added to a vertex binding offset.
MAX_VERTEX_ATTRIB_RELATIVE_OFFSET = GL_MAX_VERTEX_ATTRIB_RELATIVE_OFFSET,
// returns a single integer value containing the maximum number of vertex buffers that may be bound.
MAX_VERTEX_ATTRIB_BINDINGS = GL_MAX_VERTEX_ATTRIB_BINDINGS,
// When used with non-indexed variants of glGet (such as glGetIntegerv), returns four values: the x and y window coordinates of the viewport, followed by its width and height. Initially the x and y window coordinates are both set to 0, and the width and height are set to the width and height of the window into which the GL will do its rendering. See glViewport. When used with indexed variants of glGet (such as glGetIntegeri_v), returns four values: the x and y window coordinates of the indexed viewport, followed by its width and height. Initially the x and y window coordinates are both set to 0, and the width and height are set to the width and height of the window into which the GL will do its rendering. See glViewportIndexedf.
VIEWPORT = GL_VIEWPORT,
// returns two values, the minimum and maximum viewport bounds range. The minimum range should be at least [-32768, 32767].
VIEWPORT_BOUNDS_RANGE = GL_VIEWPORT_BOUNDS_RANGE,
// returns one value, the implementation dependent specifc vertex of a primitive that is used to select the viewport index. If the value returned is equivalent to GL_PROVOKING_VERTEX, then the vertex selection follows the convention specified by glProvokingVertex. If the value returned is equivalent to GL_FIRST_VERTEX_CONVENTION, then the selection is always taken from the first vertex in the primitive. If the value returned is equivalent to GL_LAST_VERTEX_CONVENTION, then the selection is always taken from the last vertex in the primitive. If the value returned is equivalent to GL_UNDEFINED_VERTEX, then the selection is not guaranteed to be taken from any specific vertex in the primitive.
VIEWPORT_INDEX_PROVOKING_VERTEX = GL_VIEWPORT_INDEX_PROVOKING_VERTEX,
// returns a single value, the number of bits of sub-pixel precision which the GL uses to interpret the floating point viewport bounds. The minimum value is 0.
VIEWPORT_SUBPIXEL_BITS = GL_VIEWPORT_SUBPIXEL_BITS,
