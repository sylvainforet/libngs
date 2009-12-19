/**
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>

#include "fastq.h"


static char *fastq_parser_name = NULL;

int fastq_qual_char_2_int[] = 
{
  0 , /* 0   */
  0 , /* 1   */
  0 , /* 2   */
  0 , /* 3   */
  0 , /* 4   */
  0 , /* 5   */
  0 , /* 6   */
  0 , /* 7   */
  0 , /* 8   */
  0 , /* 9   */
  0 , /* 10  */
  0 , /* 11  */
  0 , /* 12  */
  0 , /* 13  */
  0 , /* 14  */
  0 , /* 15  */
  0 , /* 16  */
  0 , /* 17  */
  0 , /* 18  */
  0 , /* 19  */
  0 , /* 20  */
  0 , /* 21  */
  0 , /* 22  */
  0 , /* 23  */
  0 , /* 24  */
  0 , /* 25  */
  0 , /* 26  */
  0 , /* 27  */
  0 , /* 28  */
  0 , /* 29  */
  0 , /* 30  */
  0 , /* 31  */
  0 , /* 32  */
  0 , /* 33  */
  0 , /* 34  */
  0 , /* 35  */
  0 , /* 36  */
  0 , /* 37  */
  0 , /* 38  */
  0 , /* 39  */
  0 , /* 40  */
  0 , /* 41  */
  0 , /* 42  */
  0 , /* 43  */
  0 , /* 44  */
  0 , /* 45  */
  0 , /* 46  */
  0 , /* 47  */
  0 , /* 48  */
  0 , /* 49  */
  0 , /* 50  */
  0 , /* 51  */
  0 , /* 52  */
  0 , /* 53  */
  0 , /* 54  */
  0 , /* 55  */
  0 , /* 56  */
  0 , /* 57  */
  0 , /* 58  */
  0 , /* 59  */
  0 , /* 60  */
  0 , /* 61  */
  0 , /* 62  */
  0 , /* 63  */
  0 , /* 64  -- @ */
  1 , /* 65  -- A */
  2 , /* 66  -- B */
  3 , /* 67  -- C */
  4 , /* 68  -- D */
  5 , /* 69  -- E */
  6 , /* 70  -- F */
  7 , /* 71  -- G */
  8 , /* 72  -- H */
  9 , /* 73  -- I */
  10, /* 74  -- J */
  11, /* 75  -- K */
  12, /* 76  -- L */
  13, /* 77  -- M */
  14, /* 78  -- N */
  15, /* 79  -- O */
  16, /* 80  -- P */
  17, /* 81  -- Q */
  18, /* 82  -- R */
  19, /* 83  -- S */
  20, /* 84  -- T */
  21, /* 85  -- U */
  22, /* 86  -- V */
  23, /* 87  -- W */
  24, /* 88  -- X */
  25, /* 89  -- Y */
  26, /* 90  -- Z */
  27, /* 91  -- [ */
  28, /* 92  -- \ */
  29, /* 93  -- ] */
  30, /* 94  -- ^ */
  31, /* 95  -- _ */
  32, /* 96  -- ` */
  33, /* 97  -- a */
  34, /* 98  -- b */
  35, /* 99  -- c */
  36, /* 100 -- d */
  37, /* 101 -- e */
  38, /* 102 -- f */
  39, /* 103 -- g */
  0 , /* 104 */
  0 , /* 105 */
  0 , /* 106 */
  0 , /* 107 */
  0 , /* 108 */
  0 , /* 109 */
  0 , /* 110 */
  0 , /* 111 */
  0 , /* 112 */
  0 , /* 113 */
  0 , /* 114 */
  0 , /* 115 */
  0 , /* 116 */
  0 , /* 117 */
  0 , /* 118 */
  0 , /* 119 */
  0 , /* 120 */
  0 , /* 121 */
  0 , /* 122 */
  0 , /* 123 */
  0 , /* 124 */
  0 , /* 125 */
  0 , /* 126 */
  0   /* 127 */
};

char *fastq_qual_char_2_string[] = 
{
  " 0 ", /* 0   */
  " 0 ", /* 1   */
  " 0 ", /* 2   */
  " 0 ", /* 3   */
  " 0 ", /* 4   */
  " 0 ", /* 5   */
  " 0 ", /* 6   */
  " 0 ", /* 7   */
  " 0 ", /* 8   */
  " 0 ", /* 9   */
  " 0 ", /* 10  */
  " 0 ", /* 11  */
  " 0 ", /* 12  */
  " 0 ", /* 13  */
  " 0 ", /* 14  */
  " 0 ", /* 15  */
  " 0 ", /* 16  */
  " 0 ", /* 17  */
  " 0 ", /* 18  */
  " 0 ", /* 19  */
  " 0 ", /* 20  */
  " 0 ", /* 21  */
  " 0 ", /* 22  */
  " 0 ", /* 23  */
  " 0 ", /* 24  */
  " 0 ", /* 25  */
  " 0 ", /* 26  */
  " 0 ", /* 27  */
  " 0 ", /* 28  */
  " 0 ", /* 29  */
  " 0 ", /* 30  */
  " 0 ", /* 31  */
  " 0 ", /* 32  */
  " 0 ", /* 33  */
  " 0 ", /* 34  */
  " 0 ", /* 35  */
  " 0 ", /* 36  */
  " 0 ", /* 37  */
  " 0 ", /* 38  */
  " 0 ", /* 39  */
  " 0 ", /* 40  */
  " 0 ", /* 41  */
  " 0 ", /* 42  */
  " 0 ", /* 43  */
  " 0 ", /* 44  */
  " 0 ", /* 45  */
  " 0 ", /* 46  */
  " 0 ", /* 47  */
  " 0 ", /* 48  */
  " 0 ", /* 49  */
  " 0 ", /* 50  */
  " 0 ", /* 51  */
  " 0 ", /* 52  */
  " 0 ", /* 53  */
  " 0 ", /* 54  */
  " 0 ", /* 55  */
  " 0 ", /* 56  */
  " 0 ", /* 57  */
  " 0 ", /* 58  */
  " 0 ", /* 59  */
  " 0 ", /* 60  */
  " 0 ", /* 61  */
  " 0 ", /* 62  */
  " 0 ", /* 63  */
  " 0 ", /* 64  -- @ */
  " 1 ", /* 65  -- A */
  " 2 ", /* 66  -- B */
  " 3 ", /* 67  -- C */
  " 4 ", /* 68  -- D */
  " 5 ", /* 69  -- E */
  " 6 ", /* 70  -- F */
  " 7 ", /* 71  -- G */
  " 8 ", /* 72  -- H */
  " 9 ", /* 73  -- I */
  " 10", /* 74  -- J */
  " 11", /* 75  -- K */
  " 12", /* 76  -- L */
  " 13", /* 77  -- M */
  " 14", /* 78  -- N */
  " 15", /* 79  -- O */
  " 16", /* 80  -- P */
  " 17", /* 81  -- Q */
  " 18", /* 82  -- R */
  " 19", /* 83  -- S */
  " 20", /* 84  -- T */
  " 21", /* 85  -- U */
  " 22", /* 86  -- V */
  " 23", /* 87  -- W */
  " 24", /* 88  -- X */
  " 25", /* 89  -- Y */
  " 26", /* 90  -- Z */
  " 27", /* 91  -- [ */
  " 28", /* 92  -- \ */
  " 29", /* 93  -- ] */
  " 30", /* 94  -- ^ */
  " 31", /* 95  -- _ */
  " 32", /* 96  -- ` */
  " 33", /* 97  -- a */
  " 34", /* 98  -- b */
  " 35", /* 99  -- c */
  " 36", /* 100 -- d */
  " 37", /* 101 -- e */
  " 38", /* 102 -- f */
  " 39", /* 103 -- g */
  " 0 ", /* 104 */
  " 0 ", /* 105 */
  " 0 ", /* 106 */
  " 0 ", /* 107 */
  " 0 ", /* 108 */
  " 0 ", /* 109 */
  " 0 ", /* 110 */
  " 0 ", /* 111 */
  " 0 ", /* 112 */
  " 0 ", /* 113 */
  " 0 ", /* 114 */
  " 0 ", /* 115 */
  " 0 ", /* 116 */
  " 0 ", /* 117 */
  " 0 ", /* 118 */
  " 0 ", /* 119 */
  " 0 ", /* 120 */
  " 0 ", /* 121 */
  " 0 ", /* 122 */
  " 0 ", /* 123 */
  " 0 ", /* 124 */
  " 0 ", /* 125 */
  " 0 ", /* 126 */
  " 0 "  /* 127 */
};

static void iter_fastq_simple (char         *path,
                               FastqIterFunc func,
                               void         *data,
                               GError      **error);

static void iter_fastq_ugly   (char         *path,
                               FastqIterFunc func,
                               void         *data,
                               GError      **error);

FastqSeq*
fastq_seq_new (void)
{
  FastqSeq *fastq;

  fastq = g_slice_new0 (FastqSeq);

  return fastq;
}

void
fastq_seq_free (FastqSeq *fastq)
{
  if (fastq)
    {
      if (fastq->name)
        g_free (fastq->name);
      if (fastq->seq)
        g_free (fastq->seq);
      if (fastq->qual)
        g_free (fastq->qual);
      g_slice_free (FastqSeq, fastq);
    }
}

void
iter_fastq (char         *path,
            FastqIterFunc func,
            void         *data,
            GError      **error)
{
  if (fastq_parser_name == NULL || g_strcmp0 (fastq_parser_name, "simple") == 0)
    iter_fastq_simple (path, func, data, error);
  else if (g_strcmp0 (fastq_parser_name, "ugly") == 0)
    iter_fastq_ugly (path, func, data, error);
  else
    {
      g_printerr ("[ERROR] Unknown fastq parser: %s\n",
                  fastq_parser_name);
      /* TODO raise an error instead of exiting */
      exit (1);
    }
}

static void
iter_fastq_simple (char         *path,
                   FastqIterFunc func,
                   void         *data,
                   GError      **error)
{
  GIOChannel *channel;
  GError     *tmp_err = NULL;
  FastqSeq   *fastq   = NULL;
  char       *line    = NULL;
  gsize       length;
  gsize       endl;
  int         ret;

  if (path[0] == '-' && path[1] == '\0')
    {
      channel = g_io_channel_unix_new (STDIN_FILENO);
      if (!channel) /* TODO raise an error here */
        return;
    }
  else
    {
      channel = g_io_channel_new_file (path, "r", &tmp_err);

      if (tmp_err != NULL)
        {
          g_propagate_error (error, tmp_err);
          if (channel)
            g_io_channel_unref (channel);
          return;
        }
    }

  while (G_IO_STATUS_NORMAL == g_io_channel_read_line (channel, &line, &length, &endl, &tmp_err))
    {
      if (*line == '@')
        {
          /* Sequence header */
          line[endl]  = '\0';
          fastq       = fastq_seq_new ();
          fastq->name = line;

          /* Sequence */
          g_io_channel_read_line (channel, &line, &length, &endl, &tmp_err);
          if (tmp_err)
            goto error;
          line[endl] = '\0';
          fastq->seq = line;

          /* Quality header */
          g_io_channel_read_line (channel, &line, &length, &endl, &tmp_err);
          if (tmp_err)
            goto error;
          g_free (line);

          /* Quality */
          g_io_channel_read_line (channel, &line, &length, &endl, &tmp_err);
          if (tmp_err)
            goto error;
          line[endl]  = '\0';
          fastq->qual = line;

          /* Size */
          fastq->size = MIN (strlen (fastq->seq),
                             strlen (fastq->qual));

          /* Callback */
          ret = func (fastq, data);
          fastq_seq_free (fastq);
          if (!ret)
            break;
        }
    }
  goto cleanup;

error:
  if (line)
    g_free (line);
  if (fastq)
    fastq_seq_free (fastq);
  g_propagate_error (error, tmp_err);

cleanup:
  g_io_channel_unref (channel);
}

static void
iter_fastq_ugly  (char         *path,
                   FastqIterFunc func,
                   void         *data,
                   GError      **error)
{
#define BUFFER_SIZE (1024 * 16)
#define MAX_STRING_SIZE 256
  char      buffer[BUFFER_SIZE];
  char      name[MAX_STRING_SIZE + 1];
  char      seq[MAX_STRING_SIZE + 1];
  char      qual[MAX_STRING_SIZE + 1];
  FILE     *file;
  FastqSeq  current;
  FastqSeq *fastq    = NULL;
  int       name_idx = 0;
  int       seq_idx  = 0;
  int       qual_idx = 0;
  int       status   = 0;
  int       i;

  if (path[0] == '-' && path[1] == '\0')
    file = stdin;
  else
    {
      file = fopen (path, "r");
      if (file == NULL)
        {
          g_printerr ("[ERROR] Could not open %s\n", path);
          /* TODO Should raise an error instead => have to create my own error domain ... */
          exit (1);
        }
    }

  do
    {
      const size_t bytes_read = fread (buffer, sizeof (char), BUFFER_SIZE, file);

      for (i = 0; i < bytes_read; i++)
        {
          if (status == 0 && buffer[i] == '@')
            {
              if (fastq)
                {
                  int ret;

                  fastq->size = MIN (seq_idx, qual_idx);
                  ret         = func (fastq, data);
                  if (!ret)
                    goto cleanup;
                }
              else
                {
                  fastq       = &current;
                  fastq->name = name;
                  fastq->seq  = seq;
                  fastq->qual = qual;
                }
              name_idx = 0;
              seq_idx  = 0;
              qual_idx = 0;
              status   = 1;
            }
          /* name */
          if (status == 1)
            {
              if (buffer[i] == '\n')
                {
                  name[name_idx] = '\0';
                  status = 2;
                }
              else if (name_idx == MAX_STRING_SIZE)
                {
                  if (buffer[i] == '\n')
                    {
                      name[name_idx] = '\0';
                      status = 2;
                    }
                }
              else
                {
                  name[name_idx] = buffer[i];
                  ++name_idx;
                }
            }
          /* sequence */
          else if (status == 2)
            {
              if (buffer[i] == '\n')
                {
                  seq[seq_idx] = '\0';
                  status = 3;
                }
              else if (seq_idx == MAX_STRING_SIZE)
                {
                  if (buffer[i] == '\n')
                    {
                      seq[seq_idx] = '\0';
                      status = 3;
                    }
                }
              else
                {
                  seq[seq_idx] = buffer[i];
                  ++seq_idx;
                }
            }
          /* name again */
          else if (status == 3)
            {
              if (buffer[i] == '\n')
                status = 4;
            }
          /* quality */
          else if (status == 4)
            {
              if (buffer[i] == '\n')
                {
                  qual[qual_idx] = '\0';
                  status = 0;
                }
              else if (qual_idx == MAX_STRING_SIZE)
                {
                  if (buffer[i] == '\n')
                    {
                      qual[qual_idx] = '\0';
                      status = 0;
                    }
                }
              else
                {
                  qual[qual_idx] = buffer[i];
                  ++qual_idx;
                }
            }
        }
      if (bytes_read < BUFFER_SIZE)
        {
          if (feof (file))
            break;
          if (ferror (file))
            {
              g_printerr ("[ERROR] Failed while reading %s\n", path);
              /* TODO Should raise an error instead */
              exit (1);
            }
        }
    }
  while (1);
  if ((status == 0 || status == 4) && fastq)
    {
      fastq->size = MIN (seq_idx, qual_idx);
      func (fastq, data);
    }

cleanup:
  fclose (file);

#undef BUFFER_SIZE
#undef MAX_NAME_SIZE
}

GOptionGroup*
get_fastq_option_group (void)
{
  GOptionEntry entries[] =
    {
      {"fastq_parser_name", 0, 0, G_OPTION_ARG_STRING, &fastq_parser_name, "Name of the parser", NULL},
      {NULL}
    };
  GOptionGroup *option_group;

  option_group = g_option_group_new ("fastq",
                                     "fastq parser options",
                                     "Show fastq parser options",
                                     NULL,
                                     NULL);
  g_option_group_add_entries (option_group, entries);

  return option_group;
}

/* vim:ft=c:expandtab:sw=4:ts=4:sts=4:cinoptions={.5s^-2n-2(0:
 */
