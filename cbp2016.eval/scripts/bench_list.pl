#!/usr/bin/perl -w
#*************************************************************
# (C) COPYRIGHT 2016 Samsung Electronics
#
#*************************************************************
# Benchmark Sets
# ************************************************************
# Provided Benchmarks are divided into 4 groups : Short_Mobile , Short_Server, Long_Mobile and Long_Server depending on the number of instructions


%SUITES = ();


$SUITES{'random'}      =
'SHORT_SERVER-32
SHORT_SERVER-40
SHORT_SERVER-45
SHORT_SERVER-77
SHORT_SERVER-85
SHORT_SERVER-99';

$SUITES{'eval'}      =
'LONG_MOBILE-1
LONG_MOBILE-10
LONG_MOBILE-11
LONG_MOBILE-12
LONG_MOBILE-13
LONG_MOBILE-14
LONG_MOBILE-15
LONG_MOBILE-16
LONG_MOBILE-17
LONG_MOBILE-18
LONG_MOBILE-19
LONG_MOBILE-2
LONG_MOBILE-3
LONG_MOBILE-4
LONG_MOBILE-5
LONG_MOBILE-6
LONG_MOBILE-7
LONG_MOBILE-8
LONG_MOBILE-9
LONG_SERVER-1
LONG_SERVER-2
LONG_SERVER-3
LONG_SERVER-4
SHORT_MOBILE-1
SHORT_MOBILE-10
SHORT_MOBILE-11
SHORT_MOBILE-12
SHORT_MOBILE-13
SHORT_MOBILE-14
SHORT_MOBILE-15
SHORT_MOBILE-16
SHORT_MOBILE-17
SHORT_MOBILE-18
SHORT_MOBILE-19
SHORT_MOBILE-2
SHORT_MOBILE-20
SHORT_MOBILE-21
SHORT_MOBILE-22
SHORT_MOBILE-23
SHORT_MOBILE-24
SHORT_MOBILE-25
SHORT_MOBILE-26
SHORT_MOBILE-27
SHORT_MOBILE-28
SHORT_MOBILE-29
SHORT_MOBILE-3
SHORT_MOBILE-30
SHORT_MOBILE-31
SHORT_MOBILE-32
SHORT_MOBILE-33
SHORT_MOBILE-34
SHORT_MOBILE-35
SHORT_MOBILE-36
SHORT_MOBILE-37
SHORT_MOBILE-38
SHORT_MOBILE-39
SHORT_MOBILE-4
SHORT_MOBILE-40
SHORT_MOBILE-41
SHORT_MOBILE-42
SHORT_MOBILE-43
SHORT_MOBILE-44
SHORT_MOBILE-45
SHORT_MOBILE-46
SHORT_MOBILE-47
SHORT_MOBILE-48
SHORT_MOBILE-49
SHORT_MOBILE-5
SHORT_MOBILE-50
SHORT_MOBILE-51
SHORT_MOBILE-52
SHORT_MOBILE-53
SHORT_MOBILE-54
SHORT_MOBILE-55
SHORT_MOBILE-56
SHORT_MOBILE-57
SHORT_MOBILE-58
SHORT_MOBILE-59
SHORT_MOBILE-6
SHORT_MOBILE-60
SHORT_MOBILE-61
SHORT_MOBILE-7
SHORT_MOBILE-8
SHORT_MOBILE-9
SHORT_SERVER-1
SHORT_SERVER-10
SHORT_SERVER-100
SHORT_SERVER-101
SHORT_SERVER-102
SHORT_SERVER-103
SHORT_SERVER-104
SHORT_SERVER-105
SHORT_SERVER-106
SHORT_SERVER-107
SHORT_SERVER-108
SHORT_SERVER-109
SHORT_SERVER-11
SHORT_SERVER-110
SHORT_SERVER-111
SHORT_SERVER-112
SHORT_SERVER-113
SHORT_SERVER-114
SHORT_SERVER-115
SHORT_SERVER-116
SHORT_SERVER-117
SHORT_SERVER-118
SHORT_SERVER-119
SHORT_SERVER-12
SHORT_SERVER-120
SHORT_SERVER-121
SHORT_SERVER-122
SHORT_SERVER-123
SHORT_SERVER-124
SHORT_SERVER-125
SHORT_SERVER-126
SHORT_SERVER-127
SHORT_SERVER-128
SHORT_SERVER-129
SHORT_SERVER-13
SHORT_SERVER-130
SHORT_SERVER-131
SHORT_SERVER-132
SHORT_SERVER-133
SHORT_SERVER-134
SHORT_SERVER-135
SHORT_SERVER-136
SHORT_SERVER-137
SHORT_SERVER-138
SHORT_SERVER-139
SHORT_SERVER-14
SHORT_SERVER-15
SHORT_SERVER-16
SHORT_SERVER-17
SHORT_SERVER-18
SHORT_SERVER-19
SHORT_SERVER-2
SHORT_SERVER-20
SHORT_SERVER-21
SHORT_SERVER-22
SHORT_SERVER-23
SHORT_SERVER-24
SHORT_SERVER-25
SHORT_SERVER-26
SHORT_SERVER-27
SHORT_SERVER-28
SHORT_SERVER-29
SHORT_SERVER-3
SHORT_SERVER-30
SHORT_SERVER-31
SHORT_SERVER-32
SHORT_SERVER-33
SHORT_SERVER-34
SHORT_SERVER-35
SHORT_SERVER-36
SHORT_SERVER-37
SHORT_SERVER-38
SHORT_SERVER-39
SHORT_SERVER-4
SHORT_SERVER-40
SHORT_SERVER-41
SHORT_SERVER-42
SHORT_SERVER-43
SHORT_SERVER-44
SHORT_SERVER-45
SHORT_SERVER-46
SHORT_SERVER-47
SHORT_SERVER-48
SHORT_SERVER-49
SHORT_SERVER-5
SHORT_SERVER-50
SHORT_SERVER-51
SHORT_SERVER-52
SHORT_SERVER-53
SHORT_SERVER-54
SHORT_SERVER-55
SHORT_SERVER-56
SHORT_SERVER-57
SHORT_SERVER-58
SHORT_SERVER-59
SHORT_SERVER-6
SHORT_SERVER-60
SHORT_SERVER-61
SHORT_SERVER-62
SHORT_SERVER-63
SHORT_SERVER-64
SHORT_SERVER-65
SHORT_SERVER-66
SHORT_SERVER-67
SHORT_SERVER-68
SHORT_SERVER-69
SHORT_SERVER-7
SHORT_SERVER-70
SHORT_SERVER-71
SHORT_SERVER-72
SHORT_SERVER-73
SHORT_SERVER-74
SHORT_SERVER-75
SHORT_SERVER-76
SHORT_SERVER-77
SHORT_SERVER-78
SHORT_SERVER-79
SHORT_SERVER-8
SHORT_SERVER-80
SHORT_SERVER-81
SHORT_SERVER-82
SHORT_SERVER-83
SHORT_SERVER-84
SHORT_SERVER-85
SHORT_SERVER-86
SHORT_SERVER-87
SHORT_SERVER-88
SHORT_SERVER-89
SHORT_SERVER-9
SHORT_SERVER-90
SHORT_SERVER-91
SHORT_SERVER-92
SHORT_SERVER-93
SHORT_SERVER-94
SHORT_SERVER-95
SHORT_SERVER-96
SHORT_SERVER-97
SHORT_SERVER-98
SHORT_SERVER-99
';


$SUITES{'SHORT_SERVER'}      =
'SHORT_SERVER-100
SHORT_SERVER-101
SHORT_SERVER-102
SHORT_SERVER-103
SHORT_SERVER-104
SHORT_SERVER-105
SHORT_SERVER-106
SHORT_SERVER-107
SHORT_SERVER-108
SHORT_SERVER-109
SHORT_SERVER-10
SHORT_SERVER-110
SHORT_SERVER-111
SHORT_SERVER-112
SHORT_SERVER-113
SHORT_SERVER-114
SHORT_SERVER-115
SHORT_SERVER-116
SHORT_SERVER-117
SHORT_SERVER-118
SHORT_SERVER-119
SHORT_SERVER-11
SHORT_SERVER-120
SHORT_SERVER-121
SHORT_SERVER-122
SHORT_SERVER-123
SHORT_SERVER-124
SHORT_SERVER-125
SHORT_SERVER-126
SHORT_SERVER-127
SHORT_SERVER-128
SHORT_SERVER-129
SHORT_SERVER-12
SHORT_SERVER-130
SHORT_SERVER-131
SHORT_SERVER-132
SHORT_SERVER-133
SHORT_SERVER-134
SHORT_SERVER-135
SHORT_SERVER-136
SHORT_SERVER-137
SHORT_SERVER-138
SHORT_SERVER-139
SHORT_SERVER-13
SHORT_SERVER-140
SHORT_SERVER-141
SHORT_SERVER-142
SHORT_SERVER-143
SHORT_SERVER-144
SHORT_SERVER-145
SHORT_SERVER-146
SHORT_SERVER-147
SHORT_SERVER-148
SHORT_SERVER-149
SHORT_SERVER-14
SHORT_SERVER-150
SHORT_SERVER-151
SHORT_SERVER-152
SHORT_SERVER-153
SHORT_SERVER-154
SHORT_SERVER-155
SHORT_SERVER-156
SHORT_SERVER-157
SHORT_SERVER-158
SHORT_SERVER-159
SHORT_SERVER-15
SHORT_SERVER-160
SHORT_SERVER-161
SHORT_SERVER-162
SHORT_SERVER-163
SHORT_SERVER-164
SHORT_SERVER-165
SHORT_SERVER-166
SHORT_SERVER-167
SHORT_SERVER-168
SHORT_SERVER-169
SHORT_SERVER-16
SHORT_SERVER-170
SHORT_SERVER-171
SHORT_SERVER-172
SHORT_SERVER-173
SHORT_SERVER-174
SHORT_SERVER-175
SHORT_SERVER-176
SHORT_SERVER-177
SHORT_SERVER-178
SHORT_SERVER-179
SHORT_SERVER-17
SHORT_SERVER-180
SHORT_SERVER-181
SHORT_SERVER-182
SHORT_SERVER-183
SHORT_SERVER-184
SHORT_SERVER-185
SHORT_SERVER-186
SHORT_SERVER-187
SHORT_SERVER-188
SHORT_SERVER-189
SHORT_SERVER-18
SHORT_SERVER-190
SHORT_SERVER-191
SHORT_SERVER-192
SHORT_SERVER-193
SHORT_SERVER-194
SHORT_SERVER-195
SHORT_SERVER-196
SHORT_SERVER-197
SHORT_SERVER-198
SHORT_SERVER-199
SHORT_SERVER-19
SHORT_SERVER-1
SHORT_SERVER-200
SHORT_SERVER-201
SHORT_SERVER-202
SHORT_SERVER-203
SHORT_SERVER-204
SHORT_SERVER-205
SHORT_SERVER-206
SHORT_SERVER-207
SHORT_SERVER-208
SHORT_SERVER-209
SHORT_SERVER-20
SHORT_SERVER-210
SHORT_SERVER-211
SHORT_SERVER-212
SHORT_SERVER-213
SHORT_SERVER-214
SHORT_SERVER-215
SHORT_SERVER-216
SHORT_SERVER-217
SHORT_SERVER-218
SHORT_SERVER-219
SHORT_SERVER-21
SHORT_SERVER-220
SHORT_SERVER-221
SHORT_SERVER-222
SHORT_SERVER-223
SHORT_SERVER-224
SHORT_SERVER-225
SHORT_SERVER-226
SHORT_SERVER-227
SHORT_SERVER-228
SHORT_SERVER-229
SHORT_SERVER-22
SHORT_SERVER-230
SHORT_SERVER-231
SHORT_SERVER-232
SHORT_SERVER-233
SHORT_SERVER-234
SHORT_SERVER-235
SHORT_SERVER-236
SHORT_SERVER-237
SHORT_SERVER-238
SHORT_SERVER-239
SHORT_SERVER-23
SHORT_SERVER-240
SHORT_SERVER-241
SHORT_SERVER-242
SHORT_SERVER-243
SHORT_SERVER-244
SHORT_SERVER-245
SHORT_SERVER-246
SHORT_SERVER-247
SHORT_SERVER-248
SHORT_SERVER-249
SHORT_SERVER-24
SHORT_SERVER-250
SHORT_SERVER-251
SHORT_SERVER-252
SHORT_SERVER-253
SHORT_SERVER-254
SHORT_SERVER-255
SHORT_SERVER-256
SHORT_SERVER-257
SHORT_SERVER-258
SHORT_SERVER-259
SHORT_SERVER-25
SHORT_SERVER-260
SHORT_SERVER-261
SHORT_SERVER-262
SHORT_SERVER-263
SHORT_SERVER-264
SHORT_SERVER-265
SHORT_SERVER-266
SHORT_SERVER-267
SHORT_SERVER-268
SHORT_SERVER-269
SHORT_SERVER-26
SHORT_SERVER-270
SHORT_SERVER-271
SHORT_SERVER-272
SHORT_SERVER-273
SHORT_SERVER-274
SHORT_SERVER-275
SHORT_SERVER-276
SHORT_SERVER-277
SHORT_SERVER-278
SHORT_SERVER-279
SHORT_SERVER-27
SHORT_SERVER-280
SHORT_SERVER-281
SHORT_SERVER-282
SHORT_SERVER-283
SHORT_SERVER-284
SHORT_SERVER-285
SHORT_SERVER-286
SHORT_SERVER-287
SHORT_SERVER-288
SHORT_SERVER-289
SHORT_SERVER-28
SHORT_SERVER-290
SHORT_SERVER-291
SHORT_SERVER-292
SHORT_SERVER-293
SHORT_SERVER-29
SHORT_SERVER-2
SHORT_SERVER-30
SHORT_SERVER-31
SHORT_SERVER-32
SHORT_SERVER-33
SHORT_SERVER-34
SHORT_SERVER-35
SHORT_SERVER-36
SHORT_SERVER-37
SHORT_SERVER-38
SHORT_SERVER-39
SHORT_SERVER-3
SHORT_SERVER-40
SHORT_SERVER-41
SHORT_SERVER-42
SHORT_SERVER-43
SHORT_SERVER-44
SHORT_SERVER-45
SHORT_SERVER-46
SHORT_SERVER-47
SHORT_SERVER-48
SHORT_SERVER-49
SHORT_SERVER-4
SHORT_SERVER-50
SHORT_SERVER-51
SHORT_SERVER-52
SHORT_SERVER-53
SHORT_SERVER-54
SHORT_SERVER-55
SHORT_SERVER-56
SHORT_SERVER-57
SHORT_SERVER-58
SHORT_SERVER-59
SHORT_SERVER-5
SHORT_SERVER-60
SHORT_SERVER-61
SHORT_SERVER-62
SHORT_SERVER-63
SHORT_SERVER-64
SHORT_SERVER-65
SHORT_SERVER-66
SHORT_SERVER-67
SHORT_SERVER-68
SHORT_SERVER-69
SHORT_SERVER-6
SHORT_SERVER-70
SHORT_SERVER-71
SHORT_SERVER-72
SHORT_SERVER-73
SHORT_SERVER-74
SHORT_SERVER-75
SHORT_SERVER-76
SHORT_SERVER-77
SHORT_SERVER-78
SHORT_SERVER-79
SHORT_SERVER-7
SHORT_SERVER-80
SHORT_SERVER-81
SHORT_SERVER-82
SHORT_SERVER-83
SHORT_SERVER-84
SHORT_SERVER-85
SHORT_SERVER-86
SHORT_SERVER-87
SHORT_SERVER-88
SHORT_SERVER-89
SHORT_SERVER-8
SHORT_SERVER-90
SHORT_SERVER-91
SHORT_SERVER-92
SHORT_SERVER-93
SHORT_SERVER-94
SHORT_SERVER-95
SHORT_SERVER-96
SHORT_SERVER-97
SHORT_SERVER-98
SHORT_SERVER-99
SHORT_SERVER-9';

$SUITES{'SHORT_MOBILE'}      =
'SHORT_MOBILE-100
SHORT_MOBILE-101
SHORT_MOBILE-102
SHORT_MOBILE-103
SHORT_MOBILE-104
SHORT_MOBILE-105
SHORT_MOBILE-106
SHORT_MOBILE-107
SHORT_MOBILE-10
SHORT_MOBILE-11
SHORT_MOBILE-12
SHORT_MOBILE-13
SHORT_MOBILE-14
SHORT_MOBILE-15
SHORT_MOBILE-16
SHORT_MOBILE-17
SHORT_MOBILE-18
SHORT_MOBILE-19
SHORT_MOBILE-1
SHORT_MOBILE-20
SHORT_MOBILE-21
SHORT_MOBILE-22
SHORT_MOBILE-23
SHORT_MOBILE-24
SHORT_MOBILE-25
SHORT_MOBILE-26
SHORT_MOBILE-27
SHORT_MOBILE-28
SHORT_MOBILE-29
SHORT_MOBILE-2
SHORT_MOBILE-30
SHORT_MOBILE-31
SHORT_MOBILE-32
SHORT_MOBILE-33
SHORT_MOBILE-34
SHORT_MOBILE-35
SHORT_MOBILE-36
SHORT_MOBILE-37
SHORT_MOBILE-38
SHORT_MOBILE-39
SHORT_MOBILE-3
SHORT_MOBILE-40
SHORT_MOBILE-41
SHORT_MOBILE-42
SHORT_MOBILE-43
SHORT_MOBILE-44
SHORT_MOBILE-45
SHORT_MOBILE-46
SHORT_MOBILE-47
SHORT_MOBILE-48
SHORT_MOBILE-49
SHORT_MOBILE-4
SHORT_MOBILE-50
SHORT_MOBILE-51
SHORT_MOBILE-52
SHORT_MOBILE-53
SHORT_MOBILE-54
SHORT_MOBILE-55
SHORT_MOBILE-56
SHORT_MOBILE-57
SHORT_MOBILE-58
SHORT_MOBILE-59
SHORT_MOBILE-5
SHORT_MOBILE-60
SHORT_MOBILE-61
SHORT_MOBILE-62
SHORT_MOBILE-63
SHORT_MOBILE-64
SHORT_MOBILE-65
SHORT_MOBILE-66
SHORT_MOBILE-67
SHORT_MOBILE-68
SHORT_MOBILE-69
SHORT_MOBILE-6
SHORT_MOBILE-70
SHORT_MOBILE-71
SHORT_MOBILE-72
SHORT_MOBILE-73
SHORT_MOBILE-74
SHORT_MOBILE-75
SHORT_MOBILE-76
SHORT_MOBILE-77
SHORT_MOBILE-78
SHORT_MOBILE-79
SHORT_MOBILE-7
SHORT_MOBILE-80
SHORT_MOBILE-81
SHORT_MOBILE-82
SHORT_MOBILE-83
SHORT_MOBILE-84
SHORT_MOBILE-85
SHORT_MOBILE-86
SHORT_MOBILE-87
SHORT_MOBILE-88
SHORT_MOBILE-89
SHORT_MOBILE-8
SHORT_MOBILE-90
SHORT_MOBILE-91
SHORT_MOBILE-92
SHORT_MOBILE-93
SHORT_MOBILE-94
SHORT_MOBILE-95
SHORT_MOBILE-96
SHORT_MOBILE-97
SHORT_MOBILE-98
SHORT_MOBILE-99
SHORT_MOBILE-9';


$SUITES{'LONG_MOBILE'}      =
'LONG_MOBILE-10
LONG_MOBILE-11
LONG_MOBILE-12
LONG_MOBILE-13
LONG_MOBILE-14
LONG_MOBILE-15
LONG_MOBILE-16
LONG_MOBILE-17
LONG_MOBILE-18
LONG_MOBILE-19
LONG_MOBILE-1
LONG_MOBILE-20
LONG_MOBILE-21
LONG_MOBILE-22
LONG_MOBILE-23
LONG_MOBILE-24
LONG_MOBILE-25
LONG_MOBILE-26
LONG_MOBILE-27
LONG_MOBILE-28
LONG_MOBILE-29
LONG_MOBILE-2
LONG_MOBILE-30
LONG_MOBILE-31
LONG_MOBILE-32
LONG_MOBILE-3
LONG_MOBILE-4
LONG_MOBILE-5
LONG_MOBILE-6
LONG_MOBILE-7
LONG_MOBILE-8
LONG_MOBILE-9';

$SUITES{'LONG_SERVER'}      =
'LONG_SERVER-1
LONG_SERVER-2
LONG_SERVER-3
LONG_SERVER-4
LONG_SERVER-5
LONG_SERVER-6
LONG_SERVER-7
LONG_SERVER-8';





$SUITES{'short'}      =
            "$SUITES{'SHORT_SERVER'}"."\n".
              "$SUITES{'SHORT_MOBILE'}";


$SUITES{'long'}      =
	      "$SUITES{'LONG_SERVER'}"."\n".
          "$SUITES{'LONG_MOBILE'}";


$SUITES{'all'}      =
	      "$SUITES{'LONG_SERVER'}"."\n".
          "$SUITES{'LONG_MOBILE'}"."\n".
            "$SUITES{'SHORT_SERVER'}"."\n".
              "$SUITES{'SHORT_MOBILE'}";
