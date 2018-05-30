using System;
using System.Collections.Generic;
using System.Collections.Concurrent;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Diagnostics;


namespace dllForMatlab
{
    public class InternalFill
    {

        byte r;
        byte m_allowed;
        List<int[]> m_dir = new List<int[]>();
        List<byte[,,]> m_hit_dist = new List<byte[,,]>();
        List<byte[,,]> m_hit = new List<byte[,,]>();
        byte[,,] m_num_dirs_found;
        int[] m_size;
        public InternalFill(byte[] range, byte[] misses_allowed, int[] size)
        {
            m_size = size;
            r = range[0];
            m_allowed = misses_allowed[0];


            for (int z = -1; z <= 1; ++z)
            {
                for (int y = -1; y <= 1; ++y)
                {
                    for (int x = -1; x <= 1; ++x)
                    {
                        m_dir.Add(new int[3] { x, y, z });
                    }
                }
            }


            //make dir arrs
            //m_dir = new byte[3, 16];
            for (int i = 0; i<27; ++i)
            {
                m_hit_dist.Add(new byte[size[0], size[1], size[2]]);
                m_hit.Add(new byte[size[0], size[1], size[2]]);
            }

            m_num_dirs_found = new byte[size[0], size[1], size[2]];
        }



        public byte[,,] Run(Byte[,,] padded_bone_map, Int32[] s1, Int32[] s2, Int32[] s3)//indexed into bone map
        {
            Stopwatch sw = new Stopwatch();
            Stopwatch sw_point = new Stopwatch();
            Byte One = (Byte)1;

            int dirs = m_dir.Count;
            //iterate through points
            //iterate through directions
            for (int i_tp = 0; i_tp < s1.Length; ++i_tp)
            {
                //sw_point.Start();

                var cs1 = s1[i_tp] - 1;
                var cs2 = s2[i_tp] - 1;
                var cs3 = s3[i_tp] - 1;

                var pcs1 = cs1 + r;
                var pcs2 = cs2 + r;
                var pcs3 = cs3 + r;

                for (int i_d = 0; i_d < dirs; ++i_d)
                {
                    var cur_dir = m_dir[i_d];
                    byte cur_val = m_hit_dist[i_d][cs1, cs2, cs3];
                    while (cur_val < r)
                    {
                        var b_val = padded_bone_map[pcs1 + ((cur_val & cur_dir[0])), pcs2 + ((cur_val & cur_dir[1])), pcs3 + ((cur_val & cur_dir[2]))];
                        if (b_val == 1)
                        {

                            for (int u = 0; u < cur_val; ++u)
                            {
                                //sw.Start();

                                m_hit_dist[i_d][cs1 + ((u) & cur_dir[0]), cs2 + ((u) & cur_dir[1]), cs3 + ((u) & cur_dir[2])] = (byte)(r + One);
                                //m_hit[i_d][cs1 + (u * (cur_val & cur_dir[0])), cs2 + u * ((cur_val & cur_dir[1])), cs3 + u * ((cur_val & cur_dir[2]))] = 1;

                                m_num_dirs_found[cs1 + ((u) & cur_dir[0]), cs2 + ((u) & cur_dir[1]), cs3 + ((u) & cur_dir[2])] += One;
                                //sw.Stop();
                                //var t = sw.ElapsedMilliseconds;


                            }
                            break;
                        }
                        else
                        {
                            ++cur_val;
                        }
                    }
                }
                //sw_point.Stop();
                //var t_1 = sw.ElapsedMilliseconds;
            }


            return m_num_dirs_found;

        }

        /// <summary>
        /// Uses a paralel look for the direction, which requirees seperate hit counters for each direction, therefore returns an array of 3d hit arrays
        /// </summary>
        /// <param name="padded_bone_map"></param>
        /// <param name="s1"></param>
        /// <param name="s2"></param>
        /// <param name="s3"></param>
        /// <returns></returns>
        public List<byte[,,]> RunParallel(Byte[,,] padded_bone_map, Int32[] s1, Int32[] s2, Int32[] s3)//indexed into bone map
        {
            Stopwatch sw = new Stopwatch();
            Stopwatch sw_point = new Stopwatch();
            Byte One = (Byte)1;

            int dirs = m_dir.Count;
            //iterate through points
            //iterate through directions
            Parallel.For(0, dirs - 1, i_d =>
            //for (int i_d = 0; i_d < dirs; ++i_d)
            {
                for (int i_tp = 0; i_tp < s1.Length; ++i_tp)
                {
                    //sw_point.Start();

                    var cs1 = s1[i_tp] - 1;
                    var cs2 = s2[i_tp] - 1;
                    var cs3 = s3[i_tp] - 1;

                    var pcs1 = cs1 + r;
                    var pcs2 = cs2 + r;
                    var pcs3 = cs3 + r;


                    var cur_dir = m_dir[i_d];
                    byte cur_val = m_hit_dist[i_d][cs1, cs2, cs3];
                    while (cur_val < r)
                    {
                        var b_val = padded_bone_map[pcs1 + ((cur_val & cur_dir[0])), pcs2 + ((cur_val & cur_dir[1])), pcs3 + ((cur_val & cur_dir[2]))];
                        if (b_val == 1)
                        {

                            for (int u = 0; u < cur_val; ++u)
                            {
                                //sw.Start();

                                m_hit_dist[i_d][cs1 + ((u) & cur_dir[0]), cs2 + ((u) & cur_dir[1]), cs3 + ((u) & cur_dir[2])] = (byte)(r + One);
                                m_hit[i_d][cs1 + ((u) & cur_dir[0]), cs2 + ((u) & cur_dir[1]), cs3 + ((u) & cur_dir[2])] = 1;

                                //m_num_dirs_found[cs1 + ((u) & cur_dir[0]), cs2 + ((u) & cur_dir[1]), cs3 + ((u) & cur_dir[2])] += One;
                                //sw.Stop();
                                //var t = sw.ElapsedMilliseconds;


                            }
                            break;
                        }
                        else
                        {
                            ++cur_val;
                        }
                    }
                }
                //sw_point.Stop();
                //var t_1 = sw.ElapsedMilliseconds;
            });








            return m_hit;
        }

        public int[] RunToIds(Byte[,,] padded_bone_map, Int32[] s1, Int32[] s2, Int32[] s3)//indexed into bone map
        {
            Stopwatch sw = new Stopwatch();
            Stopwatch sw_point = new Stopwatch();
            Byte One = (Byte)1;

            int dirs = m_dir.Count;
            //iterate through points
            //iterate through directions
            Parallel.For(0, dirs - 1, i_d =>
            //for (int i_d = 0; i_d < dirs; ++i_d)
            {
                for (int i_tp = 0; i_tp < s1.Length; ++i_tp)
                {
                    //sw_point.Start();

                    var cs1 = s1[i_tp] - 1;
                    var cs2 = s2[i_tp] - 1;
                    var cs3 = s3[i_tp] - 1;

                    var pcs1 = cs1 + r;
                    var pcs2 = cs2 + r;
                    var pcs3 = cs3 + r;


                    var cur_dir = m_dir[i_d];
                    byte cur_val = m_hit_dist[i_d][cs1, cs2, cs3];
                    while (cur_val < r)
                    {
                        var b_val = padded_bone_map[pcs1 + ((cur_val & cur_dir[0])), pcs2 + ((cur_val & cur_dir[1])), pcs3 + ((cur_val & cur_dir[2]))];
                        if (b_val == 1)
                        {

                            for (int u = 0; u < cur_val; ++u)
                            {
                                //sw.Start();

                                m_hit_dist[i_d][cs1 + ((u) & cur_dir[0]), cs2 + ((u) & cur_dir[1]), cs3 + ((u) & cur_dir[2])] = (byte)(r + One);
                                m_hit[i_d][cs1 + ((u) & cur_dir[0]), cs2 + ((u) & cur_dir[1]), cs3 + ((u) & cur_dir[2])] = 1;

                                //m_num_dirs_found[cs1 + ((u) & cur_dir[0]), cs2 + ((u) & cur_dir[1]), cs3 + ((u) & cur_dir[2])] += One;
                                //sw.Stop();
                                //var t = sw.ElapsedMilliseconds;


                            }
                            break;
                        }
                        else
                        {
                            ++cur_val;
                        }
                    }
                }
                //sw_point.Stop();
                //var t_1 = sw.ElapsedMilliseconds;
            });

            ConcurrentBag<int> passes = new ConcurrentBag<int>();
            Parallel.For(0, s1.Length - 1, (i) => 
            { 
            
                //sw_point.Start();

                var cs1 = s1[i] - 1;
                var cs2 = s2[i] - 1;
                var cs3 = s3[i] - 1;

                byte count = 0;
                for(int i_d = 0; i_d<dirs; ++i_d)
                {
                    count += m_hit[i_d][cs1, cs2, cs3];
                }
                if((dirs-count)<=m_allowed)
                {
                    passes.Add(i);

                }
            });





            return passes.ToArray();
        }



    }
}
