using C2CSharp.Models;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Web;

namespace C2CSharp.Service
{
    public class MSPLFFTService
    {
        /// <summary>
        /// 分裂基FFT算法
        /// </summary>
        /// <param name="realArr">转换前数据实部数组，长度必须是2的平方项</param>
        /// <param name="imagArr">转换前数据虚部数组，长度必须是2的平方项</param>
        /// <param name="isign">当isign等于-1的时，表示正向变换；当isign等于1的时，表示逆向变换</param>
        /// <returns>最终计算通过算法计算得到的结果</returns>
        public FFTResultViewModel msplfft(double[] realArr, double[] imagArr, int isign = -1)
        {
            FFTResultViewModel data = new FFTResultViewModel();//最终的返回结果载体
            if (!(realArr != null && imagArr != null && realArr.Length == imagArr.Length && realArr.Length > 0))//参数判断
            {
                data.Type = false;
                data.Message = "参数不符合规范！";
                return data;
            }

            //分裂基快速傅里叶变换
            try
            {
                ComplexViewModel xt = new ComplexViewModel();
                double es, e, a, a3, cc1, ss1, cc3, ss3, r1, r2, s1, s2, s3;
                int m, n2, k, n4, j, iss, id, i1, i2, i3, i0, n1, i;
                int n = realArr.Length;

                //判断维度是否为2的n次方，并确定n的值。
                double pow = Math.Log(n, 2);
                if (pow % 1 == 0)
                    m = (int)pow;
                else
                {
                    data.Type = false;
                    data.Message = "输入数据的长度必须是2的次方项！";
                    return data;
                }

                n2 = n << 1;
                es = -isign * 3.1415926535897931 * 2;
                int arrlength = n / 4;

                //存储计算的值
                double[] cc1arr = new double[arrlength];
                double[] ss1arr = new double[arrlength];
                double[] cc3arr = new double[arrlength];
                double[] ss3arr = new double[arrlength];
             
                a = 0.0;
                for (j = 0; j < n >> 2; j++)
                {
                    e = es / n;
                    a3 = 3 * a;
                    cc1arr[j] = cc1 = Math.Cos(a);
                    ss1arr[j] = ss1 = Math.Sin(a);
                    cc3arr[j] = cc3 = Math.Cos(a3);
                    ss3arr[j] = ss3 = Math.Sin(a3);
                    a = (j + 1) * e;
                }

                for (k = 1; k < m; k++)
                {
                    n2 = n2 >> 1;
                    n4 = n2 >> 2;
                    for (j = 0; j < n4; j++)
                    {

                        int idx = j * n / n2;
                        cc1 = cc1arr[idx];
                        ss1 = ss1arr[idx];
                        cc3 = cc3arr[idx];
                        ss3 = ss3arr[idx];
                        iss = j;
                        id = n2 << 1;
                        do
                        {
                            for (i0 = iss; i0 < n; i0 += id)
                            {
                                i1 = i0 + n4;
                                i2 = i1 + n4;
                                i3 = i2 + n4;
                                r1 = realArr[i0] - realArr[i2];
                                s1 = imagArr[i0] - imagArr[i2];
                                r2 = realArr[i1] - realArr[i3];
                                s2 = imagArr[i1] - imagArr[i3];
                                realArr[i0] += realArr[i2];
                                imagArr[i0] += imagArr[i2];
                                realArr[i1] += realArr[i3];
                                imagArr[i1] += imagArr[i3];
                                if (isign != 1)
                                {
                                    s3 = r1 - s2;
                                    r1 += s2;
                                    s2 = r2 - s1;
                                    r2 += s1;
                                }
                                else
                                {
                                    s3 = r1 + s2;
                                    r1 = r1 - s2;
                                    s2 = -r2 - s1;
                                    r2 = -r2 + s1;
                                }
                                realArr[i2] = r1 * cc1 - s2 * ss1;
                                imagArr[i2] = -s2 * cc1 - r1 * ss1;
                                realArr[i3] = s3 * cc3 + r2 * ss3;
                                imagArr[i3] = r2 * cc3 - s3 * ss3;
                            }

                            iss = (id << 1) - n2 + j;
                            id = id << 2;
                        } while (iss < n - 1);
                    }
                }

                /*   ------------ special last stage -------------------------*/
                iss = 0;
                id = 4;
                do
                {
                    for (i0 = iss; i0 < n; i0 += id)
                    {
                        i1 = i0 + 1;
                        xt.real = realArr[i0];
                        xt.imag = imagArr[i0];
                        realArr[i0] = xt.real + realArr[i1];
                        imagArr[i0] = xt.imag + imagArr[i1];
                        realArr[i1] = xt.real - realArr[i1];
                        imagArr[i1] = xt.imag - imagArr[i1];
                    }

                    iss = (id << 1) - 2;
                    id = id << 2;
                } while (iss < n - 1);
                j = 1;
                n1 = n - 1;
                for (i = 1; i <= n1; i++)
                {
                    if (i < j)
                    {
                        xt.real = realArr[j - 1];
                        xt.imag = imagArr[j - 1];
                        realArr[j - 1] = realArr[i - 1];
                        imagArr[j - 1] = imagArr[i - 1];
                        realArr[i - 1] = xt.real;
                        imagArr[i - 1] = xt.imag;
                    }
                    k = n >> 1;
                    do
                    {
                        if (k >= j)
                            break;
                        j -= k;
                        k = k >> 1;
                    } while (true);
                    j += k;
                }
                if (isign != -1)
                {
                    for (i = 0; i < n; i++)
                    {
                        realArr[i] /= (double)n;
                        imagArr[i] /= (double)(n);
                    }
                }

                data.CalcRealArr = realArr;
                data.CalcImagArr = imagArr;
                return data;
            }
            catch (Exception ex)
            {
                data.Type = false;
                data.Message = ex.Message;
                return data;
            }
        }
    }
}