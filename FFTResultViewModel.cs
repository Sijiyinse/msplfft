using System;
using System.Collections.Generic;
using System.Linq;
using System.Web;

namespace C2CSharp.Models
{
    /// <summary>
    /// FFT算法返回的结果ViewModel
    /// </summary>
    public class FFTResultViewModel : InfoViewModel
    {
        /// <summary>
        /// 返回的结果集
        /// </summary>
        public double[] CalcRealArr { set; get; }
        public double[] CalcImagArr { set; get; }
    }
}