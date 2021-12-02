#pragma once

#include <cstdio>
#include <sstream>
#include <boost/lexical_cast.hpp>

namespace kt84 {

namespace {
#ifdef _WIN32
    inline FILE* popen(const char* command, const char* mode) { return _popen(command, mode); }
    inline int pclose(FILE* file) { return _pclose(file); }
#endif
}

namespace zenity_util {
    inline int call(std::string& result, const std::string& param) {
        std::stringstream command;
        command << "zenity " << param;
        FILE* p = popen(command.str().c_str(), "r");
        const int BUF_SIZE = 4096;
        char buf[BUF_SIZE];
        buf[0] = 0;
        while (fgets(buf, BUF_SIZE, p) != nullptr) {}
        
        // remove '\n' stored at the back of the buffer
        int result_len = strlen(buf);
        if (result_len > 0) buf[result_len - 1] = 0;
        
        result = buf;
        
        return pclose(p);
    }
    
    // text entry
    inline bool entry(
        std::string& result,
        const std::string& title      = "",
        const std::string& text       = "",
        const std::string& entry_text = "",
        bool hide_text = false)
    {
        std::stringstream param;
        param << "--entry";
        
        if (!title      .empty()) param << " --title=\""       << title       << "\"";
        if (!text       .empty()) param << " --text=\""        << text        << "\"";
        if (!entry_text .empty()) param << " --entry-text=\""  << entry_text  << "\"";
        if (hide_text) param << " --hide-text";
        
        return call(result, param.str()) == 0;
    }
    
    // text entry for non-string type via boost::lexical_cast
    template <typename T>
    inline bool entryT(
        T& result,
        const std::string& title      = "",
        const std::string& text       = "",
        const T& init_value = T(),
        bool hide_text = false)
    {
        std::stringstream param;
        param << "--entry";
        
        if (!title      .empty()) param << " --title=\""       << title       << "\"";
        if (!text       .empty()) param << " --text=\""        << text        << "\"";
        if (hide_text) param << " --hide-text";
        
        std::string entry_text = boost::lexical_cast<std::string>(init_value);
        param << " --entry-text=\""  << entry_text  << "\"";
        
        std::string result_str;
        if (call(result_str, param.str()))
            return false;
        
        try {
            result = boost::lexical_cast<T>(result_str);
        } catch (const boost::bad_lexical_cast&) {
            return false;
        }
        return true;
    }
    
    // dialog box with OK button only
    enum class MsgBoxType {
        Error,
        Warning,
        Info,
    };
    inline void msgbox(
        MsgBoxType type,
        const std::string& title = "",
        const std::string& text  = "",
        bool no_wrap = false)
    {
        std::string result;
        std::stringstream param;
        param << "--" << (type == MsgBoxType::Error ? "error" : type == MsgBoxType::Warning ? "warning" : "info");
        
        if (!title.empty()) param << " --title=\"" << title << "\"";
        if (!text .empty()) param << " --text=\""  << text  << "\"";
        if (no_wrap) param << " --no-wrap";
        
        call(result, param.str());
    }
    
    // dialog box with OK/Cancel buttons
    inline bool question(
        const std::string& title        = "",
        const std::string& text         = "",
        const std::string& ok_label     = "",
        const std::string& cancel_label = "",
        bool no_wrap = false)
    {
        std::string result;
        std::stringstream param;
        param << "--question";
        
        if (!title       .empty()) param << "        --title=\"" << title        << "\"";
        if (!text        .empty()) param << "         --text=\"" << text         << "\"";
        if (!ok_label    .empty()) param << "     --ok-label=\"" << ok_label     << "\"";
        if (!cancel_label.empty()) param << " --cancel-label=\"" << cancel_label << "\"";
        if (no_wrap) param << " --no-wrap";
        
        return call(result, param.str()) == 0;
    }
    
    // file selection
    inline bool file_selection_load(
        std::string& result,
        const std::string& title       = "",
        const std::string& filename    = "",
        const std::string& file_filter = "",
        bool directory = false,
        bool multiple = false,
        char separator = 0)
    {
        std::stringstream param;
        param << "--file-selection";
        
        if (!title      .empty()) param << "       --title=\"" << title       << "\"";
        if (!filename   .empty()) param << "    --filename=\"" << filename    << "\"";
        if (!file_filter.empty()) param << " --file-filter=\"" << file_filter << "\"";
        if (directory) param << " --directory";
        if (multiple ) param << "  --multiple";
        if (separator) param << " --separator=" << separator;
        
        return call(result, param.str()) == 0;
    }
    
    inline bool file_selection_save(
        std::string& result,
        const std::string& title       = "",
        const std::string& filename    = "",
        const std::string& file_filter = "",
        bool confirm_overwrite = true ,
        bool directory         = false,
        bool multiple          = false,
        char separator = 0)
    {
        std::stringstream param;
        param << "--file-selection --save";
        
        if (!title      .empty()) param << "       --title=\"" << title       << "\"";
        if (!filename   .empty()) param << "    --filename=\"" << filename    << "\"";
        if (!file_filter.empty()) param << " --file-filter=\"" << file_filter << "\"";
        if (confirm_overwrite) param << " --confirm-overwrite";
        if (directory        ) param << " --directory";
        if (multiple         ) param << " --multiple";
        if (separator) param << " --separator=" << separator;
        
        return call(result, param.str()) == 0;
    }
    
    // scale bar
    inline bool scale(
        int& result,
        const std::string& title = "",
        const std::string& text  = "",
        int value       =   0,
        int min_value   =   0,
        int max_value   = 100,
        int step        =   1,
        bool hide_value = false)
    {
        std::string result_str;
        std::stringstream param;
        param << "--scale";
        if (!title.empty()) param << " --title=\"" << title << "\"";
        if (!text .empty()) param << " --text=\""  << text  << "\"";
        param << "     --value=" << value    ;
        param << " --min-value=" << min_value;
        param << " --max-value=" << max_value;
        param << "      --step=" << step     ;
        if (hide_value) param << " --hide-value";
        
        if (call(result_str, param.str()))
            return false;
        
        result = boost::lexical_cast<int>(result_str);
        return true;
    }
    
    // reporting with text box
    inline void text_info(
        std::string& result,
        const std::string& title    = "",
        const std::string& filename = "",
        bool editable = false)
    {
        /// TODO
    }
    
    //--calendar
    //--list
    //--notification
    //--progress
}

}

