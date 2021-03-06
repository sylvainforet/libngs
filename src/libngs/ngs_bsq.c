/* Copyright (C) 2010  Sylvain FORET
 *
 * This file is part of libngs.
 *
 * libngs is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *                                                                       
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *                                                                       
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/**
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>

#include "ngs_bsq.h"
#include "ngs_bsq_flex.h"
#include "ngs_utils.h"


static char *bsq_parser_name = NULL;

static void iter_bsq_simple (char         *path,
                             BsqIterFunc   func,
                             void         *data,
                             GError      **error);

char *strand_names[BSQ_STRAND_NB] =
{
  [BSQ_STRAND_W ] = "++",
  [BSQ_STRAND_WC] = "+-",
  [BSQ_STRAND_C ] = "-+",
  [BSQ_STRAND_CC] = "--"
};

char *map_flag_names[BSQ_MAP_FLAG_NB] =
{
  [BSQ_MAP_UM] = "UM",
  [BSQ_MAP_MA] = "MA",
  [BSQ_MAP_OF] = "OF",
  [BSQ_MAP_NM] = "NM",
  [BSQ_MAP_QC] = "QC"
};

BsqRecord*
bsq_record_new (void)
{
  BsqRecord *rec;

  rec = g_slice_new0 (BsqRecord);

  return rec;
}

void
bsq_record_free (BsqRecord *rec)
{
  if (rec)
    {
      if (rec->name)
        g_free (rec->name);
      if (rec->ref)
        g_free (rec->ref);
      if (rec->seq)
        g_free (rec->seq);
      if (rec->mis_info)
        g_free (rec->mis_info);
      if (rec->mC_loc)
        g_free (rec->mC_loc);
      g_slice_free (BsqRecord, rec);
    }
}

void
iter_bsq (char         *path,
          BsqIterFunc   func,
          void         *data,
          GError      **error)
{
  if (bsq_parser_name == NULL || strcmp (bsq_parser_name, "flex") == 0)
    iter_bsq_flex (path, func, data, error);
  else if (strcmp (bsq_parser_name, "simple") == 0)
    iter_bsq_simple (path, func, data, error);
  else
    {
      g_set_error (error,
                   NGS_ERROR,
                   NGS_UNKNOWN_ERROR,
                   "Unknown bsq parser: `%s'",
                   bsq_parser_name);
    }
}

static void
iter_bsq_simple (char         *path,
                 BsqIterFunc   func,
                 void         *data,
                 GError      **error)
{
  GIOChannel *channel;
  GError     *tmp_err = NULL;
  char       *line    = NULL;
  BsqRecord   record;
  gsize       length;
  gsize       endl;

  if (path[0] == '-' && path[1] == '\0')
    {
      channel = g_io_channel_unix_new (STDIN_FILENO);
      if (!channel)
        {
          g_set_error (error,
                       NGS_ERROR,
                       NGS_IO_ERROR,
                       "Could not open stdin for reading");
          return;
        }
    }
  else
    {
      channel = g_io_channel_new_file (path, "r", &tmp_err);

      if (tmp_err)
        {
          g_propagate_error (error, tmp_err);
          if (channel)
            g_io_channel_unref (channel);
          return;
        }
    }

  while (G_IO_STATUS_NORMAL == g_io_channel_read_line (channel, &line, &length, &endl, &tmp_err))
    {
      int    ret;
      char **fields;

      line[endl] = '\0';
      fields = g_strsplit (line, "\t", -1);
      if (fields[0])
        {
          record.name = fields[0];
          if (fields[1])
            {
              record.seq = fields[1];
              record.size = strlen (fields[1]);
              if (fields[2])
                {
                  switch (fields[2][0])
                    {
                      case 'U':
                          record.flag = BSQ_MAP_UM;
                          break;
                      case 'M':
                          record.flag = BSQ_MAP_MA;
                          break;
                      case 'O':
                          record.flag = BSQ_MAP_OF;
                          break;
                      case 'N':
                          record.flag = BSQ_MAP_NM;
                          break;
                      case 'Q':
                          record.flag = BSQ_MAP_QC;
                          break;
                      default:
                          break;
                    }
                  if (fields[3])
                    {
                      record.ref = fields[3];
                      if (fields[4])
                        {
                          record.loc = atol (fields[4]);
                          if (fields[5])
                            {
                              if (fields[5][0] == '+')
                                {
                                  if (fields[5][1] == '+')
                                    record.strand = BSQ_STRAND_W;
                                  else
                                    record.strand = BSQ_STRAND_WC;
                                }
                              else
                                {
                                  if (fields[5][1] == '+')
                                    record.strand = BSQ_STRAND_C;
                                  else
                                    record.strand = BSQ_STRAND_CC;
                                }
                              if (fields[6])
                                {
                                  record.n_mis = atoi (fields[6]);
                                  if (fields[7])
                                    {
                                      record.mis_info = fields[7];
                                      if (fields[8])
                                        {
                                          record.mC_loc = fields[8];
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
      ret = func (&record, data);
      g_strfreev (fields);
      g_free (line);
      line   = NULL;
      if (!ret)
        break;
      record.name     = NULL;
      record.seq      = NULL;
      record.ref      = NULL;
      record.mis_info = NULL;
      record.mC_loc   = NULL;
      record.size     = 0;
    }
  if (line)
    g_free (line);
  if (tmp_err)
    g_propagate_error (error, tmp_err);
  g_io_channel_unref (channel);
}

GOptionGroup*
get_bsq_option_group (void)
{
  GOptionEntry entries[] =
    {
      {"bsq_parser_name", 0, 0, G_OPTION_ARG_STRING, &bsq_parser_name, "Name of the bsq parser", NULL},
      {NULL}
    };
  GOptionGroup *option_group;

  option_group = g_option_group_new ("bsq",
                                     "bsq parser options",
                                     "Show bsq parser options",
                                     NULL,
                                     NULL);
  g_option_group_add_entries (option_group, entries);

  return option_group;
}

/* vim:ft=c:expandtab:sw=4:ts=4:sts=4:cinoptions={.5s^-2n-2(0:
 */
