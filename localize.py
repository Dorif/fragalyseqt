# This file is part of FragalyseQt.
#
# FragalyseQt is free software: you can redistribute it and/or modify it under
# the terms of the GNU Affero General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version.
#
# FragalyseQt is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License
# for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with FragalyseQt. If not, see <https://www.gnu.org/licenses/>.

def localizefq(iface):
    from os import name, getenv
    if name == 'posix':
        lang = getenv('LANG')
    elif name == 'nt':
        from locale import windows_locale
        from ctypes import windll
        lang = windows_locale[windll.kernel32.GetUserDefaultUILanguage()]
    else:
        lang = "en"
# Well, some systems like FreeBSD may have C or POSIX locale.
# Or no locale set, like OpenBSD.
    if lang is None or ("en" or "C.UTF-8" or "POSIX") in lang:
        iface['ch_inact_msg'] = "Inactive channel"
        iface['aboutbtn'] = "About"
        iface['infoboxtxt'] = "FragalyseQt version 0.4.1, codename\"Jeffreys"\
                              "\".\n\nThis program version supports "\
                              "analysis of up to 8 different fluorescense "\
                              "channels simultaneously, selective channel "\
                              "hiding, non-Latin run names and can correctly"\
                              " handle damaged files, exporting peaks "\
                              "locations, areas, sizes, FWHM's and channel "\
                              "names in CSV for any *.FSA/*.HID files and "\
                              "export in CSV internal analysis data for "\
                              "*.FSA/*.HID files generated by ABI 3500 and "\
                              "SeqStudio series equipment.\n\nPeaks areas are"\
                              " calculated assuming they are Gaussian peaks"\
                              ".\n\nLicensed under GNU AGPL version 3.\n\nIf "\
                              "you wish to contactauthor for any reason - "\
                              "please write at dorif11@gmail.com"
        iface['openfragmentfile'] = "Open file"
        iface['exportinternal'] = "Export Internal Analysis Data"
        iface['csvexport'] = "Export in CSV"
        iface['unsuppeq'] = "Unsupported equipment!"
        iface['unsuppeqmsg'] = "Internal analysis data could be exported only"\
                               " from FSA files generated by ABI 3500 and "\
                               "SeqStudio family sequencers!"
        iface['dmgdfile'] = "Damaged file!"
        iface['nodatamsg'] = "There are no data in file!"
        iface['hidechannel'] = "Hide "
        iface['minph'] = "Select minimal peak height:"
        iface['minpw'] = "Select minimal peak width:"
        iface['minpp'] = "Select minimal peak prominence:"
        iface['minww'] = "Select minimal window width:"
        iface['savecsv'] = "Save CSV"
        iface['bcd'] = "Baseline correction and denoising."
        iface['wrongsizing'] = "Wrong ladder or sizing method! Please, try "\
                               "another ones!"
    elif "ru" in lang:
        iface['ch_inact_msg'] = "Неактивный канал"
        iface['aboutbtn'] = "О программе"
        iface['infoboxtxt'] = "FragalyseQt версия 0.4.1, кодовое имя "\
                              "\"Джеффрис\".\n\nЭта версия программы "\
                              "поддерживает одновременный анализ до 8 каналов"\
                              " флуоресценции, избирательное сокрытие каналов"\
                              ", не-латинские имена запусков, может правильно"\
                              " работать с повреждёнными файлами, "\
                              "поддерживает экспорт имён каналов, положений, "\
                              "площадей, размеров и FWHM пиков в формат CSV "\
                              "для любого файла *.FSA/*.HID и экспорт данных"\
                              " внутреннего анализа данных для файлов "\
                              "*.FSA/*.HID, полученных на инструментах ABI "\
                              "3500 и SeqStudio.\n\nПлощади пиков "\
                              "рассчитываются исходя из того, что пики "\
                              "являются гауссовыми.\n\nЛицензировано под "\
                              "GNU AGPL версии 3.\n\nЕсли по какой-либо "\
                              "причине вы хотите проконтактировать с автором"\
                              " - пишите, пожалуйста, на dorif11@gmail.com"
        iface['openfragmentfile'] = "Открыть файл"
        iface['exportinternal'] = "Экспорт Данных Внутреннего Анализа"
        iface['csvexport'] = "Экспорт в CSV"
        iface['unsuppeq'] = "Неподдерживаемое оборудование!"
        iface['unsuppeqmsg'] = "Данные внутреннего анализа могут быть "\
                               "экспортированы только из файлов FSA, "\
                               "полученных на ABI 3500 или семействе "\
                               "SeqStudio!"
        iface['dmgdfile'] = "Повреждённый файл!"
        iface['nodatamsg'] = "В файле отсутствуют данные!"
        iface['hidechannel'] = "Скрыть "
        iface['minph'] = "Выберите минимальную высоту пика:"
        iface['minpw'] = "Выберите минимальную ширину пика:"
        iface['minpp'] = "Выберите минимальное возвышение пика:"
        iface['minww'] = "Выберите минимальную ширину окна:"
        iface['savecsv'] = "Сохранить CSV"
        iface['bcd'] = "Коррекция базовой линии и шумоподавление."
        iface['wrongsizing'] = "Ошибочный размерный маркер или метод "\
                               "определения размера! Пожалуйста, "\
                               "попробуйте другие!"
    elif "ua" in lang:
        iface['ch_inact_msg'] = "Неактивний канал"
        iface['aboutbtn'] = "О програмі"
        iface['infoboxtxt'] = "FragalyseQt версія 0.4.1, кодове ім'я "\
                              "\"Джефріс\".\n\n Версія програми підтримує "\
                              "одночасний аналіз до 8 каналів флуоресценции, "\
                              "вибіркове приховування каналів, не-латинські "\
                              "імена запусків, може правильно працювати з "\
                              "пошкодженими файлами, підтримує експорт імен "\
                              "каналів, положень, площадей, розміру і FWHM "\
                              "піков у формат CSV для будь-якого файла"\
                              " *.FSA/*.HID та експорт даних внутрішнього "\
                              "аналізу даних для файлів *.FSA/*.HID, "\
                              "отриманих на інструментах ABI 3500 и "\
                              "SeqStudio.\n\nПлощаді піков розраховуються "\
                              "виходячи з того, що піки є гауссовими.\n\n"\
                              "Ліцензовано під GNU AGPL версії 3.\n\nЯкщо з "\
                              "будь-якої причини ви хочете проконтактувати з "\
                              "автором - пишіть, будь ласка, на "\
                              "dorif11@gmail.com"
        iface['openfragmentfile'] = "Відкрити файл"
        iface['exportinternal'] = "Експорт Даних Внутренього Аналіза"
        iface['csvexport'] = "Експорт в CSV"
        iface['unsuppeq'] = "Непідтримуване оболаднення!"
        iface['unsuppeqmsg'] = "Дані внутренього аналіза можуть бути "\
                               "експортовані тільки із файлів FSA, отриманих"\
                               " на ABI 3500 або сімействі SeqStudio!"
        iface['dmgdfile'] = "Пошкоджений файл!"
        iface['nodatamsg'] = "У файлі відсутні дані!"
        iface['hidechannel'] = "Сховати "
        iface['minph'] = "Виберіть мінімальну висоту піку:"
        iface['minpw'] = "Виберіть мінімальну ширину піку:"
        iface['minpp'] = "Виберіть мінімальне підвищення піку:"
        iface['minww'] = "Виберіть мінімальну ширину вікна:"
        iface['savecsv'] = "Зберегти CSV"
        iface['bcd'] = "Корекція базової лінії та шумозаглушення."
        iface['wrongsizing'] = "Помилковий розмірний маркер чи метод "\
                               "визначення розміру! Будь ласка, спробуйте"\
                               " інші!"
    elif "ro" in lang:
        iface['ch_inact_msg'] = "Canal inactiv"
        iface['aboutbtn'] = "Despre program"
        iface['infoboxtxt'] = "FragalyseQt versiune 0.4.1, codename "\
                              "\"Jeffreys\".\n\nAceastă versiune permite "\
                              "analiza simultană a pînă la 8 canale de "\
                              "fluorescență, ascunderea selectivă a canalelor"\
                              ", nume non-latine a runurilor, poate funcționa"\
                              " corect cu fișiere deteriorate, permite "\
                              "exportarea numelor de canale, locațiile, "\
                              "suprafețele, dimensiunile și FWHMurile a "\
                              "vîrfurilor în formatul CSV, pentru orice "\
                              "fișierul de tip *.FSA/*.HID și exportarea în "\
                              "formatul CSV a datelor analizei interne "\
                              "pentru fișiere de tip *.FSA/*.HID generate de"\
                              " echipamentele de model "\
                              "ABI 3500 și SeqStudio.\n\nSuprafețele "\
                              "vîrfurilor, sunt calculate presupunînd că "\
                              "vîrfurile sunt gausiene.\n\nLicențiat sub GNU "\
                              "AGPL versiunea 3.\n\nDacă doriți să contactați"\
                              " autorul, pentru orice întrebare - scrieți, vă"\
                              " rog, la adresa dorif11@gmail.com"
        iface['openfragmentfile'] = "Deschideți fișierul"
        iface['exportinternal'] = "Exportul Datelor Analizei Interne"
        iface['csvexport'] = "Exportul în CSV"
        iface['unsuppeq'] = "Echipament nesusținut!"
        iface['unsuppeqmsg'] = "Datele analizei interne pot fi exportate "\
                               "numai pentru fișierile FSA generate de "\
                               "echipamentele de model ABI 3500 și SeqStudio!"
        iface['dmgdfile'] = "Fișier deteriorat!"
        iface['nodatamsg'] = "In fișier lipsesc datele!"
        iface['hidechannel'] = "Ascunde "
        iface['minph'] = "Selectaţi înălţimea minimală a vîrfului:"
        iface['minpw'] = "Selectaţi lăţimea minimală a vîrfului:"
        iface['minpp'] = "Selectaţi proeminenţa minimală a vîrfului:"
        iface['minww'] = "Selectaţi lăţimea minimală a ferestrei:"
        iface['savecsv'] = "Salvează CSV"
        iface['bcd'] = "Corectarea liniei de bază și supresia zgomotului."
        iface['wrongsizing'] = "Marker de dimensiuni sau metodă de "\
                               "dimensionare incorectă! Vă rugăm să"\
                               " încercați altele!"
    elif "fr" in lang:
        iface['ch_inact_msg'] = "Canal inactif"
        iface['aboutbtn'] = "Au sujet de"
        iface['infoboxtxt'] = "FragalyseQt version 0.4.1, alias "\
                              "\"Jeffreys\".\n\nCette version du programme "\
                              "prend en charge l'analyse simultanée de "\
                              "jusqu'à 8 canaux de fluorescence différents, "\
                              "le masquage sélectif des canaux, les noms "\
                              "d'exécution non latins et peut gérer "\
                              "correctement les Fichiers endommagés, exporter"\
                              " les emplacements des pics, zones, tailles, "\
                              "FWHM et noms de canal au format CSV pour tous "\
                              "les Fichiers *.FSA/*.HID et exportation au "\
                              "format CSV des données d'analyse interne pour "\
                              "les Fichiers *.FSA/*.HID générés par les "\
                              "équipements des séries ABI 3500 et SeqStudio."\
                              "\n\nLes zones de pics sont calculées en "\
                              "supposant qu'il s'agit de pics gaussiens.\n\n"\
                              "Licence sous GNU AGPL version 3.\n\nSi vous "\
                              "souhaitez contacter l'auteur pour une raison "\
                              "quelconque, veuillez écrire à dorif11@gmail.com"
        iface['openfragmentfile'] = "Ouvrer le Fichier"
        iface['exportinternal'] = "Exporter des données d'analyse interne"
        iface['csvexport'] = "Exporter en CSV"
        iface['unsuppeq'] = "L'équipement non supportés!"
        iface['unsuppeqmsg'] = "Les données d'analyse internes ne pouvaient "\
                               "être exportées qu'à partir de Fichiers FSA "\
                               "générés par les séquenceurs de la famille "\
                               "ABI 3500 et SeqStudio!"
        iface['dmgdfile'] = "Fichier endommagé!"
        iface['nodatamsg'] = "Il n'y a pas de données dans le Fichier!"
        iface['hidechannel'] = "Masquer"
        iface['minph'] = "Sélectionnez la hauteur de pic minimale:"
        iface['minpw'] = "Sélectionnez la largeur de pic minimale:"
        iface['minpp'] = "Sélectionnez la proéminence de pic minimale:"
        iface['minww'] = "Sélectionnez la largeur de la fenêtre minimale:"
        iface['savecsv'] = "Enregistrer CSV"
        iface['bcd'] = "Correction de la ligne de base et suppression du"\
                       " bruit."
        iface['wrongsizing'] = "Marqueur de dimensionnement ou méthode de "\
                               "dimensionnement incorrect ! S'il vous plaît"\
                               ", essayez-en d'autres !"
    elif "bg" in lang:
        iface['ch_inact_msg'] = "Неактивен канал"
        iface['aboutbtn'] = "Относно програмата"
        iface['infoboxtxt'] = "FragalyseQt версия 0.4.1, кодово име \"Джефрис"\
                              "\".\n\nТази версия на програмата поддържа "\
                              "едновременен анализ до 8 флуоресцентни канала,"\
                              " селективно скриване на канали, не-латински "\
                              "имена на изпълнения, може да работи правилно с"\
                              " повредени файлове, поддържа експорт на имена "\
                              "на канали, позиции, площите, размери и FWHM "\
                              "пикове в CSV формат за всеки *.FSA/*.HID файл "\
                              "и експортиране на вътрешни данни за анализ на "\
                              "данни за *.FSA/*.HID файлове, получени с ABI "\
                              "3500 и инструменти SeqStudio.\n\nПлощите на "\
                              "пиковете се изчисляват, като се приема, че "\
                              "пиковете са гаусови.\n\nЛицензирано под GNU "\
                              "AGPL версия 3.\n\nАко по някаква причина "\
                              "искате да се свържете с автора, моля, "\
                              "изпратете имейл на dorif11@gmail.com"
        iface['openfragmentfile'] = "Отвори файла"
        iface['exportinternal'] = "Експорт на Данни от Вътрешен Анализ"
        iface['csvexport'] = "Експорт в CSV"
        iface['unsuppeq'] = "Неподдържан оборудване!"
        iface['unsuppeqmsg'] = "Данните за вътрешен анализ могат да бъдат "\
                               "експортирани само от файлове FSA, придобити "\
                               "на ABI 3500 или фамилията SeqStudio!"
        iface['dmgdfile'] = "Повреден файл!"
        iface['nodatamsg'] = "Във файла няма данни!"
        iface['hidechannel'] = "Крия "
        iface['minph'] = "Изберете минимална височина на върха:"
        iface['minpw'] = "Изберете минимална ширина на върха:"
        iface['minpp'] = "Изберете минимална надморска височина:"
        iface['minww'] = "Изберете минимална ширина на прозореца:"
        iface['savecsv'] = "Запазване на CSV"
        iface['bcd'] = "Корекция на базовата линия и шумопотискане."
        iface['wrongsizing'] = "Неправилен маркер за оразмеряване или метод"\
                               " за оразмеряване! Моля, опитайте други!"
