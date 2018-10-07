#include <jni.h>
#include <string>
#include<math.h>
#include <cpu-features.h>
#include <vector>
#define PI 3.1415926535897
std::vector<std::vector<jdouble>> k_value;

std::vector<jdouble> kernel_values(jdouble s, jint r,jint y)
{
    jint k=0,i;
    jdouble sum = 0.0;
    jint kSize = abs(2 * r + 1);
    std::vector<jdouble> kernel_value;
    //kernel values
    for( i=0;i<kSize;i++)
    {
        k=i-r;
        jdouble gc = pow((1.0 / (2.0 * PI * s * s)), 0.5);
        jdouble gec = exp((-pow(k, 2)) /(2.0 * pow(s, 2)));
        jdouble val = gc*gec;
        kernel_value.push_back(val);
        sum+=val;
    }
    //Normalisation
    for(i=0; i<kSize; i++){
        jdouble v = kernel_value[i]/sum;
        kernel_value[i] = v;
    }
    return kernel_value;
}


void kernel_compute(jint a0 , jint a1, jint a2, jint a3,jint h, jint r, jdouble sigmanear, jdouble sigmafar) {
    jint i,e;
    jdouble s;
    std::vector<jdouble> k;
    //map of kernel values for first region
    for(i=0;i<a0;i++)
    {
        k_value.push_back(kernel_values(sigmafar, r, 0));
    }
    //map of kernel values for second region
    for (i = a0; i < a1; i++) {
        s = abs(sigmafar * (jdouble)(a1-i)/(a0-a1));
        k_value.push_back(kernel_values(s, r, i));
    }
    for(i=a1;i<a2;i++)
    {
        std::vector<jdouble> k1(3,1);
        k_value.push_back(k1);
    }
    //map of kernel values for fourth region
    for (i = a2; i < a3; i++) {
        s = sigmanear * (jdouble)((i - a2) / (a3 - a2));
        k_value.push_back(kernel_values(s, r, i));
    }
    //map of kernel values for fifth region
    for(i=a3;i<h;i++)
    {
        k_value.push_back(kernel_values(sigmanear, r, i));
    }
    int a=i;
}

std::vector<jdouble> intermediate_matrix(std::vector<jdouble> &i_matrix,std::vector<jdouble> &input, jint a0, jint a1,jint a2,jint a3, jint w, jint h, jint r)
{
    jint i,j,kernelIndex,pixelIndex,k,e=0;
    std::vector<jdouble> kernel;
    jint size=w*h;
    jdouble value;
    for(i=0;i<h;i++)
    {
        if(i<a1||i>a2) {
            // kernel selection
            kernel = k_value[i];

            // intermediate pixel computation
            for (j = 0; j < w; j++) {
                value = 0.0;
                for (k = 0; k < (2 * r) + 1; k++) {
                    pixelIndex = (i * w) + j -k;
                    if (pixelIndex >= 0 && pixelIndex < size) {
                        value += (input[pixelIndex] * kernel[kernelIndex]);

                    }
                }
                i_matrix[(i * w) + j] = value;
            }
        }
    }

return i_matrix;
}

std::vector<jdouble> output_matrix(std::vector<jdouble> &o_mat,std::vector<jdouble> &i_mat,jint a0, jint a1,jint a2,jint a3, jint w, jint h, jint r)
{
    jint kernelIndex,pixelIndex,e=0;
    std::vector<jdouble> kernel;
    for(int i=0;i<h;i++)
    {
        if(i<a1||i>a2)
        {
        // kernel computation
        kernel = k_value[i];
        // output pixel computation
            for (int j = 0; j < w; j++)
            {
                jdouble value = 0.0;
                for (jint k = 0; k < (2 * r) + 1; k++)
                {
                    kernelIndex = k;
                    pixelIndex = (i - kernelIndex) * w + j;
                    if (pixelIndex >= 0 && pixelIndex < i_mat.size())
                    {
                            value += (kernel[kernelIndex] * i_mat[pixelIndex]);

                    }
                }
                o_mat[(i * w) + j] = value;
            }
        }
    }
    return o_mat;
}

std::vector<jdouble> Build_Blur(std::vector<jdouble> &channel_input, jint a0, jint a1,jint a2,jint a3,jint w,jint h,jfloat sigma_far,jfloat sigma_near)
{
    jint radius= 10;
    std::vector<jdouble> output,int_matrix;
    output = channel_input;
    int_matrix = channel_input;


    kernel_compute(a0,a1,a2, a3,h,radius,(jdouble)sigma_near,(jdouble)sigma_far);
    int_matrix = intermediate_matrix(int_matrix,channel_input,a0,a1,a2,a3,w,h,radius);
    output = output_matrix(output,int_matrix,a0,a1,a2,a3,w,h,radius);
    channel_input.clear();
    return output;
}

extern "C"
JNIEXPORT jint JNICALL
Java_edu_asu_ame_meteor_speedytiltshift2018_SpeedyTiltShift_tiltshiftcppnative(JNIEnv *env,
                                                                               jobject instance,
                                                                               jintArray inputPixels_,
                                                                               jintArray outputPixels_,
                                                                               jint width,
                                                                               jint height,
                                                                               jfloat sigma_far,
                                                                               jfloat sigma_near,
                                                                               jint a0, jint a1,
                                                                               jint a2, jint a3) {
    // Declarations
    jint *pixels = env->GetIntArrayElements(inputPixels_, NULL);
    jint *outputPixels = env->GetIntArrayElements(outputPixels_, NULL);
    jint size = width*height;
    jint a=0xff;
    int i;
    std::vector<jdouble> channelR;
    std::vector<jdouble> channelG;
    std::vector<jdouble> channelB;
    std::vector<jdouble> output_channelR;
    std::vector<jdouble> output_channelG;
    std::vector<jdouble> output_channelB;

    // create a copy of input pixels and three channels
    for (int j=0;j<size;j++)
    {
        jint B = pixels[j] & 0xff;
        jint G = (pixels[j]>>8) & 0xff;
        jint R = (pixels[j]>>16)&0xff;
        channelR.push_back((jdouble)R);
        channelG.push_back((jdouble)G);
        channelB.push_back((jdouble)B);
    }

    //Build Blur
    output_channelR = Build_Blur(channelR,a0,a1,a2,a3,width,height,sigma_far,sigma_near);
    output_channelG = Build_Blur(channelG,a0,a1,a2,a3,width,height,sigma_far,sigma_near);
    output_channelB = Build_Blur(channelB,a0,a1,a2,a3,width,height,sigma_far,sigma_near);

    //Build final Output
    for(i=0;i<size;i++) {
        jint aVal = (a & 0xff) << 24;
        jint rVal = ((jint)output_channelR[i] & 0xff) << 16;
        jint gVal = ((jint)output_channelG[i] & 0xff) << 8;
        jint bVal = ((jint)output_channelB[i] & 0xff);
        jint val= aVal|rVal|gVal|bVal;
        outputPixels[i] = val ;
    }


    env->ReleaseIntArrayElements(inputPixels_, pixels, 0);
    env->ReleaseIntArrayElements(outputPixels_, outputPixels, 0);
     return 0;
}




extern "C"
JNIEXPORT jint JNICALL
Java_edu_asu_ame_meteor_speedytiltshift2018_SpeedyTiltShift_tiltshiftneonnative(JNIEnv *env,
                                                                               jobject instance,
                                                                               jintArray inputPixels_,
                                                                               jintArray outputPixels_,
                                                                               jint width,
                                                                               jint height,
                                                                               jfloat sigma_far,
                                                                               jfloat sigma_near,
                                                                               jint a0, jint a1,
                                                                               jint a2, jint a3) {
    jint *pixels = env->GetIntArrayElements(inputPixels_, NULL);
    jint *outputPixels = env->GetIntArrayElements(outputPixels_, NULL);

    for (int j=0;j<height;j++){
        for (int i=0;i<width;i++) {
            int B = pixels[j*width+i]%0xff;
            int G = (pixels[j*width+i]>>8)%0xff;
            int R = (pixels[j*width+i]>>16)%0xff;
            int A = 0xff;
            R=0;
            int color = (A & 0xff) << 24 | (R & 0xff) << 16 | (G & 0xff) << 8 | (B & 0xff);

            outputPixels[j*width+i]=color;
        }
    }

    env->ReleaseIntArrayElements(inputPixels_, pixels, 0);
    env->ReleaseIntArrayElements(outputPixels_, outputPixels, 0);
    return 0;
}