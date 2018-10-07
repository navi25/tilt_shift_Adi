#include <jni.h>
#include <string>
#include<math.h>
#include <cpu-features.h>
#include<unordered_map>
#include <vector>

#define PI 3.1415926535897
std::unordered_map<jdouble,std::vector<jdouble>> map;

inline jdouble compute_sigma(jint low_bound, jint high_bound, jdouble s, jint y, jint blur_type)
{
    switch(blur_type)
    {
        default:
            break;
        case 1:
            return s;
        case 2:
            s = abs(s * ((high_bound - y) / (low_bound - high_bound)));
            break;
        case 3:
            s = s * ((y-low_bound)/(high_bound-low_bound));
            break;
    }
    return s;
}

std::vector<jdouble> kernel_compute(jint r ,jdouble s)
{
    if(map.find(s) != map.end())
    {
        return map[s];
    }
    jint k=0;
    jdouble sum = 0.0;
    jint kSize = abs(2 * r + 1);
   // jdouble * kernel_value = new jdouble((kSize));
     std::vector<jdouble> kernel_value;
    for(int i=0;i<kSize;i++)
    {
        k=i-r;
        jdouble gc = pow((1.0 / (2.0 * PI * s * s)), 0.5);
        jdouble gec = exp((-pow(k, 2)) /(2.0 * pow(s, 2)));
        jdouble val = gc*gec;
        kernel_value.push_back(val);
        sum+=val;
    }

    //Normalisation
    for(int i=0; i<kSize; i++){
        jdouble v = kernel_value[i]/sum;
        kernel_value[i] = v;
    }
    map[s] = kernel_value;

    return kernel_value;

}

inline std::vector<jdouble> intermediate_matrix(std::vector<jdouble> &input, std::vector<jdouble> &kernel_value,jint l_bound, jint h_bound, jdouble sigma,jint width, jint h,jint radius)
{
        int k,kernelIndex,pixelIndex;
        int size = width*h;
        //jdouble* int_matrix = new jdouble[size];
        std::vector<jdouble> int_matrix(size,0);
        for(int i=l_bound;i<h_bound;i++)
        {

        for(int j=0;j<width;j++)
        {

            if (sigma > 0.6)
            {
                jdouble value = 0.0;

                for (k = 0; k < (2 * radius) + 1; k++)
                {
                    //kernelIndex = k-radius;
                    //int pixelIndex = (i*width)+j - radius;
                    kernelIndex = k;
                    pixelIndex = (i*width)+j +kernelIndex;
                    if(pixelIndex>=0 && pixelIndex<size)
                    {
                        if (kernelIndex >= 0)
                        {
                            value +=( input[pixelIndex] * kernel_value[kernelIndex]);

                        }
                    }

                }

                int_matrix[(i*width)+j]=value;
            }

        }
    }
    return int_matrix;

}
inline std::vector<jdouble> output_matrix(std::vector<jdouble> &output,std::vector<jdouble> &int_matrix,std::vector<jdouble> &kernel_value,jint l_bound,jint h_bound,jint width,jint h,jint radius)
{
    jint kernelIndex,pixelIndex;
    jint size =width*h;
    for(int i=l_bound;i<h_bound;i++)
    {
        for(int j=0;j<width;j++)
        {
//            output[j]= 0;
            jdouble value = 0.0;
            for(jint k=0;k<(2*radius)+1;k++)
            {
                //kernelIndex = k-radius;
                //jint pixelIndex = (j*width) +i-radius;
                kernelIndex = k;
                pixelIndex = (i-kernelIndex)*width+j;
                if(pixelIndex >= 0 && pixelIndex<int_matrix.size())
                {
                    if (kernelIndex >= 0)
                    {
                        value += ( kernel_value[kernelIndex] * int_matrix[pixelIndex]);

                    }
                }
            }
            output[(i*width)+j]=value;
        }
    }
    //delete []int_matrix;
    //int_matrix.clear();
    return output;
}





std::vector<jdouble> Gaussian_Blur(std::vector<jdouble> &output,std::vector<jdouble> &region_input,jint l_bound,jint h_bound,jint blur_type,jdouble sigma,jint w,jint h) {
    jint radius = 10;

    if (blur_type == 1) {
        std::vector<jdouble> kernel_value, int_matrix;
        sigma = compute_sigma(l_bound, h_bound, sigma, 0, 1);
        kernel_value = kernel_compute(radius, sigma);
        int_matrix = intermediate_matrix(region_input, kernel_value, l_bound, h_bound, sigma, w, h,
                                         radius);
        output = output_matrix(output, int_matrix, kernel_value, l_bound, h_bound, w, h, radius);
    }
    if (blur_type == 2)
    {
        for(int i=l_bound; i<h_bound-1; i++){
            std::vector<jdouble> kernel_value, int_matrix;
            sigma = compute_sigma(l_bound, h_bound, sigma, i, 2);
            kernel_value = kernel_compute(radius, sigma);
            int_matrix = intermediate_matrix(region_input, kernel_value, i, i+1, sigma, w, h,
                                             radius);
            output = output_matrix(output, int_matrix, kernel_value, i, i+1, w, h, radius);

        }

    }
    if (blur_type == 3)
    {
        for(int i=l_bound; i<h_bound-1; i++){
            std::vector<jdouble> kernel_value, int_matrix;
            sigma = compute_sigma(l_bound, h_bound, sigma, i, 3);
            kernel_value = kernel_compute(radius, sigma);
            int_matrix = intermediate_matrix(region_input, kernel_value, i, i+1, sigma, w, h,
                                             radius);
            output = output_matrix(output, int_matrix, kernel_value, i, i+1, w, h, radius);

        }

    }
    return output;

}


std::vector<jdouble> Build_Blur(std::vector<jdouble> &channel_input, jint a0, jint a1,jint a2,jint a3,jint w,jint h,jfloat sigma_far,jfloat sigma_near)
{
    jint size = w*h;
    std::vector<jdouble> output;
    //for (int j=0;j<size;j++)
    //{
         //output[j] = channel_input[j];
      //  output.push_back(channel_input[j]);
    //}
    output = channel_input;
    // First Region
    //output = Gaussian_Blur(output,channel_input,0,a0,1,(jdouble)sigma_far,w,h);
    // second region
     output =  Gaussian_Blur(output,channel_input,a0,a1,2,(jdouble)sigma_far,w,h);
    //fourth region
    //output =  Gaussian_Blur(output,channel_input,a2,a3,3,(jdouble)sigma_near,w,h);
    // fifth region
    //output =  Gaussian_Blur(output,channel_input,a3,h,1,(jdouble)sigma_near,w,h);
            //delete []channel_input;
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
    //vector<jdouble>channelR, *channelG, *channelB,*output_channelR, *output_channelG, *output_channelB;
    std::vector<jdouble> channelR;
    std::vector<jdouble> channelG;
    std::vector<jdouble> channelB;
    std::vector<jdouble> output_channelR;
    std::vector<jdouble> output_channelG;
    std::vector<jdouble> output_channelB;
    //channelR = new jdouble[size];
    //channelG = new jdouble[size];
    //channelB = new jdouble[size];


    // create a copy of input pixels and three channels
    for (int j=0;j<size;j++)
    {
        int i = pixels[j];
        jint B = pixels[j] & 0xff;
        jint G = (pixels[j]>>8) & 0xff;
        jint R = (pixels[j]>>16)&0xff;
        channelR.push_back((jdouble)R);
        channelG.push_back((jdouble)G);
        channelB.push_back((jdouble)B);
        //channelR[j] = (jdouble)R;
        //channelG[j] = (jdouble)G;
        //channelB[j] = (jdouble)B;
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
    int k=i;
    //map.clear();
    //delete []output_channelR;
    //output_channelR.clear();
    //delete []output_channelG;
    //output_channelG.clear();
    //delete []output_channelB;
    //output_channelB.clear();

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