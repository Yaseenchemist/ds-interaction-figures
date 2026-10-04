# ds-interaction-figures

أتمتة كاملة لرسم تفاعلات الليكند مع البروتين (3D + 2D) بدقة نشر، باستخدام **Discovery Studio Visualizer 2025** وبايثون.
لكل ملف PDB تحصل على:

- **صورة 3D**: الليكند ball-and-stick (C سماوي، O أحمر، N أزرق، S أصفر، Cl أخضر)، والأحماض المتفاعلة فقط بخطوط رفيعة (الكربون أسود)، مع أسماء الأحماض وأرقامها والمسافات بخط واضح وخطوط تفاعل ملونة.
- **4 زوايا إضافية** لنفس المركب لاختيار الأوضح.
- **مخطط 2D** بنفس ألوان الـ 3D، مع مفتاح "Interactions" والمسافات (600 dpi).

## الفكرة
DS يرسم **الهندسة فقط** بدقة 3200×2400 بدون أي نص (لأن نص DS يتشوه عند الدقة العالية). ثم `finish_images.py` يتعلم موقع كاميرا DS من صور معايرة ويكتب الأسماء والمسافات فوق الصورة بخط حاد.

## الهيكل
```
ds_automation/
├── ds_batch_interactions.pl   ← يعمل داخل Discovery Studio
├── finish_images.py           ← يعمل في الترمنال
└── pdb/                       ← ضع هنا كل ملفات *.pdb
```
تُنشأ تلقائياً: `out/` (خام + جداول) و `out/final/` (الصور النهائية) و `out/final/views/` (الزوايا).

## المتطلبات
- Discovery Studio Visualizer 2025 (لينكس).
- Python 3 مع: `python3 -m pip install cairosvg pillow numpy`

## الإعداد (مرة واحدة)
1. انسخ المستودع إلى `~/Documents/ds_automation/` (أو غيّر `my $BASE` في أعلى ملف `.pl` و `--folder` في بايثون).
2. ضع ملفات PDB في `pdb/`. اسم الملف = اسم الصور الناتجة، استخدم أسماء بسيطة بدون مسافات (`lig1.pdb`).
   الليكند يُكتشف تلقائياً (LIG, UNL, `*0` ...) والمعادن (Zn وغيره) القريبة منه تظهر.

## التشغيل
1. افتح Discovery Studio Visualizer.
2. **File ▸ Open** واختر `ds_batch_interactions.pl` (ليس ملف pdb). السكربت يفتح كل ملفات `pdb/` واحداً بعد الآخر ويغلقها. لا تلمس DS أثناء العمل.
3. راقب `out/batch_log.txt`، وينتهي بسطر `ALL DONE`.
4. في الترمنال:
   ```bash
   python3 ~/Documents/ds_automation/finish_images.py
   ```
5. الصور النهائية في `out/final/`:

| الملف | المحتوى |
|---|---|
| `<name>_3D.png` | صورة 3D الرئيسية |
| `views/<name>_view1..4.png` | زوايا بديلة |
| `<name>_2D.png` | مخطط 2D (600 dpi) |

> كل تشغيل يعيد معالجة كل ما في `pdb/` ويكتب فوق الصور القديمة. انسخ `out/final/` قبل دفعة جديدة إن أردت الاحتفاظ بها.

## إعدادات شائعة
في أعلى `ds_batch_interactions.pl`:

| المتغير | الوظيفة |
|---|---|
| `@LIGAND_CARBON = (0,190,190)` | لون كربون الليكند |
| `@RESIDUE_CARBON = (0,0,0)` | لون كربون الأحماض |
| `$LINE_WIDTH = 2` | سمك خطوط الأحماض |
| `$SAVE_VIEWS = 1` | حفظ الزوايا الإضافية (`0` لإيقافها) |
| `@VIEWS` / `@VIEW_ROT` | الزوايا / تدوير الصورة الرئيسية |
| `$IMG_W, $IMG_H` | حجم الصورة (3200×2400) |
| `$SHOW_DS_LINES = 0` | `1` = إبقاء خطوط تفاعل DS (ثم شغّل بايثون مع `--ds-lines`) |

خيارات `finish_images.py`:
```bash
python3 finish_images.py --only lig1            # مركب واحد
python3 finish_images.py --dec 2                # مسافتان بمنزلتين عشريتين
python3 finish_images.py --res-font 0.02 --dist-font 0.017   # تكبير الخط
python3 finish_images.py --res-fmt cap          # Asp1046 بدل ASP 1046
python3 finish_images.py --carbon 0,190,190     # يجب أن يطابق لون الـ .pl
python3 finish_images.py --keep-vdw --keep-h    # إبقاء van der Waals والهيدروجين في 2D
```

## حل المشاكل
- **`calibration failed` بدون أسماء**: المعايرة لم تنجح لهذا المركب. افتح `out/calib/<name>_v0_m0.png`، يجب أن ترى 6 كرات ملونة مختلفة. إذا لم تظهر، أرسلها للمطور.
- **`no ligand found`**: لا يوجد ليكند HETATM. تأكد أن الليكند HETATM وليس ATOM.
- **`out/` فارغ**: أنت فتحت ملف pdb بدل `.pl`، أو السكربت توقف، فراجع `batch_log.txt`.
- **خط مفقود (Luxi Sans)**: لا يؤثر على النتيجة، فالخط يُستبدل بـ DejaVu Sans.
- أي خطوة فاشلة داخل DS تظهر في اللوغ بكلمة `FAILED` ولا توقف بقية الدفعة.

## ملاحظة
الخطوات (ألوان، ترتيب، تسمية) مضبوطة على DS 2025. بعض أوامر الـ API (تدوير الزوايا) تعمل حسب النسخة؛ الفشل يظهر في اللوغ فقط.
